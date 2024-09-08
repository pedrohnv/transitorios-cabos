#= Funções para ler uma netlist, montar e simular um circuito no domínio da
frequência com resultados no domínio do tempo obtido por transformada de Laplace.
Os passos no tempo e frequência são uniformes (constantes).

As primeiras linhas da netlist definem os parâmetros do tempo na simulação:
passo no tempo `dt` e número de passos no tempo `nt`. Por causa da transformada
inversa de Laplace, aproximadamente de 20% a 25% dos passos finais no tempo
serão perdidos por causa de lixo numério proporcional a `exp(t)`.

A netlist tem o seguinte estilo:
    dt valor
    nt valor
    NOME1 VALOR NÓ_1 NÓ_2 ... NÓ_N
    NOME2 VALOR NÓ_1 NÓ_2 ... NÓ_N

As primeiras letras de NOME identificam o tipo do elemento, que podem ser:
    V:fonte_de_tensão  tensão  nó_positivo  nó_negativo
    I:fonte_de_corrente  corrente  nó_positivo(entrada)  nó_negativo(saída)
    R:resistor  resistência  nó_positivo  nó_negativo
    L:indutor  indutância  nó_positivo  nó_negativo
    C:capacitor  capacitância  nó_positivo  nó_negativo
    Z:RLC  resistência,indutância,capacitância  nó_positivo  nó_negativo
    M:modelo_externo  arquivo  nó_1  ...  nó_n

Para todos os elementos, VALOR pode ser uma constante (no tempo) ou nome de um arquivo.
Os nós são strings na qual "0" é reservado para simbolizar o nó Terra.

No caso de VALOR ser dado por um arquivo, ele pode ser dos formatos CSV ou HDF5.
É esperado que ele tenha estampas de tempos ("t") frequências angulares
complexas em rad/s ("s"), além do valor ("v", independente do tipo de
elemento). Se não houver um valor correspondente ao tempo na simulação, então ele
é obtido por interpolação linear:
    `y = (y0 * (t1 - t) + y1 * (t - t0)) / (t1 - t0)`

O elemento Z sintetiza dentro de si um RLC série. Útil pra evitar ter que criar
nós intermediários. Os valores devem estar separados por vírgula, e não espaço!

Sobre o elemento modelo externo M, é esperado que ele seja, necessariamente,
um arquivo HDF5 com frequências angulares `s` e matrizes de admitância nodal
correspondentes `yn`. Os nós de M correspondem, em ordem, à posição na matriz.

Fontes de tensão aumentam o tamanho e piora o condicionamento numérico do sistema
linear, então é preferível utilizar um equivalente Norton, se possível, onde uma
fonte de tensão `V` em série com uma resistência `R` é substituído por uma fonte
de corrente `I=V/R` em paralelo com uma resistência `R`.

É possível ter uma linha de comentário estritamente a partir da terceira linha,
desde que o primeiro caracter dela seja '#' ou '%'.
=#
using CSV, DataFrames
using HDF5
using Interpolations

const NO_TERRA = "0";

include("laplace.jl")

"""Lê a netlist e faz alguns testes de sanidade. Retorna um dicionário dos
nome dos nós, fontes de tensão e a posição deles na matriz de admitância nodal,
e os parâmetros da simulação `dt` e `nt`."""
function preprocessar_netlist(netlist_file::String)
    local netlist
    open(netlist_file, "r") do file
        netlist = readlines(file)
    end
    # Ler parâmetros da simulação: `dt` e `nt`
    linha = netlist[1]
    elem = split(linha)
    if elem[1] != "dt"
        msg = "Esperava-se que a linha `$(linha)` tivesse a forma: `dt 1e-9`"
        throw(error(msg))
    end
    dt = tryparse(Float64, elem[2])
    if isnothing(dt)
        msg = "Valor inválido em `$(linha)`"
        throw(error(msg))
    end
    
    linha = netlist[2]
    elem = split(linha)
    if elem[1] != "nt"
        msg = "Esperava-se que a linha `$(linha)` tivesse a forma: `nt 1024`"
        throw(error(msg))
    end
    nt = tryparse(Int, elem[2])
    if isnothing(nt)
        msg = "Valor inválido em `$(linha)`"
        throw(error(msg))
    end
    
    # Extrair todos os nós e contar as fontes de tensão
    lista_nos = Vector{String}([])
    lista_fonte_tensao = Vector{String}([])
    linhas_ignorar = Vector{Int}([])
    for (e, elem) in enumerate(netlist)
        if elem == "" || elem == '\n' || elem[1] == '#' || elem[1] == '%' || elem[1:2] == "dt" || elem[1:2] == "nt"
            push!(linhas_ignorar, e)
            continue  # pular linhas vazias ou que são comentário
        end
        elemento = string.(split(elem))
        nome = elemento[1]
        tipo = nome[1:1]
        nodes = elemento[3:end]
        if tipo ∈ ["V", "I", "R", "L", "C", "Z", "M"]
            if tipo == "V"
                if nome ∈ lista_fonte_tensao
                    msg = "Nome repetido para fonte de tensão: $(nome)"
                    throw(error(msg))
                end
                push!(lista_fonte_tensao, "i_" * nome)
            end
            if tipo != "M"
                if nodes[1] == nodes[2]
                    msg = "Elemento $(nome) com nós repetidos (curto-circuitados)"
                    throw(error(msg))
                elseif length(nodes) > 2
                    msg = "Elemento $(nome) com mais de dois nós"
                    throw(error(msg))
                end
            end
            for no in nodes
                if no ∉ lista_nos
                    push!(lista_nos, no)
                end
            end
        else
            msg = "Elemento não identificado em `$(s)`"
            throw(error(msg))
        end
    end
    deleteat!(netlist, linhas_ignorar)
    if NO_TERRA ∉ lista_nos
        msg = "Nenhum elemento conectado ao nó Terra $(NO_TERRA)"
        @warn msg
    end
    # Remover o nó Terra
    d = findfirst([no == NO_TERRA for no in lista_nos])
    if !isnothing(d)
        splice!(lista_nos, d)
    end
    dict_potenciais = merge(Dict(NO_TERRA => 0), Dict([(lista_nos[i], i) for i = 1:length(lista_nos)]))
    n = length(dict_potenciais) - 1
    dict_correntes = Dict([(lista_fonte_tensao[i], n + i) for i = 1:length(lista_fonte_tensao)])
    dict_variaveis = merge(dict_potenciais, dict_correntes)
    @assert allunique(values(dict_variaveis))
    return netlist, dict_variaveis, dt, nt
end


abstract type ElementoAbstrato end


"""Estrutura para representar uma admitância genérica. Possui um vetor com
o índice dos nós positivo (node[1]) e negativo (nodes[2]) e um vetor com os valores
para cada frequência de interesse para evitar ler um arquivo múltiplas vezes."""
struct MatrizAdmitancia <: ElementoAbstrato
    nome::String
    valores::Array{ComplexF64, 3}
    nodes::Vector{Int}
end


"""Nó positivo (nodes[1]) é o que sai corrente da fonte."""
struct FonteCorrente <: ElementoAbstrato
    nome::String
    valores::Vector{ComplexF64}
    nodes::Vector{Int}
end


"""Nó positivo (nodes[1]) é o que sai corrente da fonte."""
struct FonteTensao <: ElementoAbstrato
    nome::String
    valores::Vector{ComplexF64}
    nodes::Vector{Int}
    corrente::Int  # posição (índice) da variável de corrente da fonte
end


"""Extrair o valor do elemento de um CSV."""
function valor_csv(file, dt, nt)
    tmax = dt * (nt - 1)
    df = CSV.read(file, DataFrame)
    local s, t, v, val
    if "s" in names(df)
        s = parse.(ComplexF64, df.s)
        val = parse.(ComplexF64, df.v)
        temp, freq_s = laplace(ones(nt), tmax, nt)
        #TODO interpolar s ao invés de dar erro e dar @warn quando interpolação é usada
        if !all(isapprox.(freq_s, s))
            msg = "As frequências angulares `s` no arquivo $(file) não condizem com as da simulação"
            throw(error(msg))
        end
    elseif "t" in names(df)
        t = df.t
        v = df.v
        if t[end] < tmax
            @warn "não há valor definido para `v(tmax)` em $(file). Interpolando com `v(tmax)=0`"
            push!(t, tmax)
            push!(v, 0.0)
        end
        li = linear_interpolation(t, v)
        y = li.(range(0, length=nt, step=dt))
        val, s = laplace(y, tmax, nt)
    else
        msg = "$(file) não possui uma coluna s ou t"
        throw(error(msg))
    end
    return val, s
end


"""Extrair o valor do elemento de um HDF5."""
function valor_h5(file, dt, nt)
    tmax = dt * (nt - 1)
    local s, t, v, val
    h5open(file, "r") do h5
        v = read(h5, "v")
        if "s" in keys(h5)
            s = read(h5, "s")
            #TODO interpolar s ao invés de dar erro e dar @warn quando interpolação é usada
            temp, freq_s = laplace(ones(nt), tmax, nt)
            if !all(isapprox.(freq_s, s))
                msg = "As frequências angulares `s` no arquivo $(file) não condizem com as da simulação"
                throw(error(msg))
            end
            val = v
        elseif "t" in keys(h5)
            t = read(h5, "t")
            if t[end] < tmax
                @warn "não há valor definido para `v(tmax)` em $(file). Interpolando com `v(tmax)=0`"
                push!(t, tmax)
                push!(v, 0.0)
            end
            li = linear_interpolation(t, v)
            y = li.(range(0, length=nt, step=dt))
            #TODO @warn quando interpolação é usada
            val, s = laplace(y, tmax, nt)
        else
            msg = "$(file) não possui uma coluna s ou t"
            throw(error(msg))
        end
    end
    return val, s
end


"""Interpreta a string de um elemento da netlist."""
function interpretar_elemento(
    linha::String,
    dict_variaveis::Dict{String, Int},
    freq_s::Vector{ComplexF64},
    dt::Float64,
    nt::Int
)
    tmax = dt * (nt - 1)
    elemento = split(linha)
    nome = elemento[1]
    tipo = nome[1:1]
    valor = elemento[2]
    nodes = elemento[3:end]
    nodes_int = [dict_variaveis[n] for n in nodes]
    nf = length(freq_s)
    local val
    if isfile(valor)
        if valor[end-3:end] == ".csv"
            val, s = valor_csv(valor, dt, nt)
        elseif valor[end-2:end] == ".h5"
            val, s = valor_h5(valor, dt, nt)
        else
            msg = "$(valor) não é um arquivo válido. Esperava-se .csv ou .h5"
            throw(error(msg))
        end
    elseif tipo == "Z"
        z = tryparse.(Float64, split(valor, ","))
        if any(isnothing.(z))
            msg = "$(valor) não é uma tripla `R,L,C` nem um arquivo válido"
            throw(error(msg))
        end
        R, L, C = z
        zc = (C > 0.0) ? 1.0 ./ (freq_s .* C) : 0.0
        val = @. R + freq_s * L + zc
    else
        v = tryparse(Float64, valor)
        if isnothing(v)
            msg = "$(valor) não é um número nem um arquivo válido"
            throw(error(msg))
        end
        if tipo == "V" || tipo == "I"
            val, s = laplace(fill(v, nt), tmax, nt)
        else
            if v <= 0.0
                msg = "Valor de $(nome) deve ser positivo."
                throw(error(msg))
            end
            val = fill(v, nf)
        end
    end

    function equivalenteMatrizAdmitancia(v)
        i1, i2 = nodes_int
        if i1 == i2
            msg = "Elemento $(nome) curto-circuitado"
            throw(error(msg))
        end
        # reduzir matriz se nó aterrado
        yn = Array{ComplexF64}(undef, 1, 1, nf)
        if i1 == 0
            nodes_int = [i2]
            yn[:] .= v[:]
        elseif i2 == 0
            nodes_int = [i1]
            yn[:] .= v[:]
        else
            yn = Array{ComplexF64}(undef, 2, 2, nf)
            yn[1,1,:] .= v
            yn[2,1,:] .= -v
            yn[1,2,:] .= -v
            yn[2,2,:] .= v
        end
        return MatrizAdmitancia(nome, yn, nodes_int)
    end
    
    if tipo == "R" || tipo == "Z"
        return equivalenteMatrizAdmitancia(1.0 ./ val)
    elseif tipo == "L"
        return equivalenteMatrizAdmitancia(1.0 ./ (freq_s .* val))
    elseif tipo == "C"
        return equivalenteMatrizAdmitancia(freq_s .* val)
    elseif tipo == "M"
        if length(nodes_int) != size(val)[1] == size(val)[2]
            msg = "O número de nós não condiz com o tamanho da matriz para o elemento $(nome)"
            throw(error(msg))
        end
        # reduzir matriz se nó aterrado
        local i0 = findfirst([i == 0 for i in nodes_int])
        while !isnothing(i0)
            deleteat!(nodes_int, i0)
            ids = collect(1 : (length(nodes_int) + 1))
            deleteat!(ids, i0)
            val = val[ids, ids, :]
            i0 = findfirst([i == 0 for i in nodes_int])
        end
        return MatrizAdmitancia(nome, val, nodes_int)
    elseif tipo == "I"
        return FonteCorrente(nome, val, nodes_int)
    elseif tipo == "V"
        return FonteTensao(nome, val, nodes_int, dict_variaveis["i_" * nome])
    else
        msg = "tipo $(tipo) não reconhecido para elemento $(nome)"
        throw(error(msg))
    end
end


"""Modifica as matrizes `yn(s)` e vetores de excitação `ib(s)`."""
function montar_sistema!(
    yn::Array{ComplexF64, 3},
    ib::Array{ComplexF64, 2},
    elementos::Vector{ElementoAbstrato}
)
    for f = 1:size(yn)[3]
        for elem in elementos
            tipo = typeof(elem)
            if tipo <: MatrizAdmitancia
                # Aqui é assumindo que a matriz do elemento já foi reduzida em 
                # caso de algum nó estar aterrado.
                v = elem.valores[:,:,f]
                n = elem.nodes
                ynview = view(yn, n, n, f)
                ynview[:,:] .+= v
            elseif tipo <: FonteCorrente
                v = elem.valores[f]
                i1, i2 = elem.nodes
                if i1 > 0
                    ib[i1, f] -= v
                end
                if i2 > 0
                    ib[i2, f] += v
                end
            elseif tipo <: FonteTensao
                v = elem.valores[f]
                i1, i2 = elem.nodes
                i3 = elem.corrente
                if i1 > 0
                    yn[i1, i3, f] += 1.0
                    yn[i3, i1, f] += 1.0
                end
                if i2 > 0
                    yn[i2, i3, f] -= 1.0
                    yn[i3, i2, f] -= 1.0
                end
                ib[i3, f] += v
            else
                msg = "tipo $(tipo) não reconhecido para elemento"
                throw(error(msg))
            end 
        end
    end
end


"""Interpreta e simula a netlist."""
function simular_netlist(netlist_file::String)
    netlist, dict_variaveis, dt, nt = preprocessar_netlist(netlist_file)
    tmax = dt * (nt - 1)
    x, freq_s = laplace(ones(nt), tmax, nt)
    nf = length(freq_s)
    elementos = [interpretar_elemento(linha, dict_variaveis, freq_s, dt, nt) for linha in netlist]
    variaveis_dict = Dict([(v, k) for (k, v) in dict_variaveis])
    nv = length(dict_variaveis) - 1  # Lembrar de desconsiderar o Terra!
    yn = zeros(ComplexF64, nv, nv, nf)
    ib = zeros(ComplexF64, nv, nf)
    montar_sistema!(yn, ib, elementos)
    for f = 1:nf
        # lembre-se: slicing em julia (e.g. M[1:3]) cria uma cópia se usada como
        # argumento de função. Use @view, se necessário.
        ynf = view(yn, :, :, f)
        ibf = view(ib, :, f)
        LAPACK.gesv!(ynf, ibf)
    end
    vout = Dict([variaveis_dict[k] => invlaplace(ib[k,:], freq_s, tmax, nt, "Hanning") for k = 1:nv])
    #vout = [invlaplace(ib[k,:], freq_s, tmax, nt, "Hanning") for k = 1:nv]
    return vout#, variaveis_dict
end

#netlist_file = "netlist_umbilical.txt"
#simular_netlist(netlist_file)