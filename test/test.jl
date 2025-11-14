include("src/cabos.jl")

# Tripolar com blindagem interna - setup in Julia

# Material properties (as given)
epsr_c = 3.31     # condutor central
epsr_b = 2.3      # blindagem, jaqueta da veia
epsr_a = 10       # última camada, envolve a armadura

rho_c = 1.7241e-8   # Ω·m, condutor central
rho_b = 2.2e-7      # Ω·m, blindagem
rho_a = 2.86e-8     # Ω·m, armadura

mur_c = 1.0     # permeabilidade magnética relativa do condutor central
mur_b = 1.0     # permeabilidade magnética relativa da blindagem
mur_a = 300.0   # permeabilidade magnética relativa da armadura

sig_s = 5.0     # condutividade do mar [S/m]
eps_s = 81.0     # permissividade relativa do mar

# Geometria do umbilical PT
rc = 1e-3 .* [9.6, 17.054, 18.054, 19.50]   # raios dos cabos SC
ra = 1e-3 .* [48, 59, 65]                   # raios da armadura

# Posições
x0, y0 = 0.0, 0.0

d1 = (rc[end] + 1e-12) / cos(deg2rad(30.0))
xa = x0
ya = y0 + d1 * sin(deg2rad(90.0))

xb = x0 + d1 * cos(deg2rad(210.0))
yb = y0 + d1 * sin(deg2rad(210.0))

xc = x0 + d1 * cos(deg2rad(330.0))
yc = y0 + d1 * sin(deg2rad(330.0))

# Núcleo e blindagem
fase_nucleo = CableComponent(0.0, rc[1], rc[2], rho_c, mur_c, 1.0, epsr_c)
fase_blindagem = CableComponent(rc[2], rc[3], rc[4], rho_b, mur_b, 1.0, epsr_b)

# Cabos coaxiais
fase_a = CoaxialCable([fase_nucleo, fase_blindagem], xa, ya)
fase_b = CoaxialCable([fase_nucleo, fase_blindagem], xb, yb)
fase_c = CoaxialCable([fase_nucleo, fase_blindagem], xc, yc)

# Armadura (PipeCable) contendo os três coaxiais
armadura = PipeCable(
    ra[1],          # radius_in
    ra[2],          # radius_ext
    ra[3],          # radius_ext_insulator
    rho_a,            # rho_c
    mur_a,            # mur_c
    1.0,              # mur_d_in
    1.0,              # mur_d_ext
    epsr_b,           # epsr_in
    epsr_a,           # epsr_ext
    cables = [fase_a, fase_b, fase_c],
    x = ra[1],          # position plausível, pode ajustar conforme necessidade
    y = 0.0,              # posição y
    name = "Pipe",
    sigma_medium = 5.0,
    epsr_medium = 81.0,
    cable_length = 1.0
)
