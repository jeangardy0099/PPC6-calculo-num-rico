import os
import time
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# CONFIGURAÇÕES GERAIS
# -----------------------------------------------------------------------------
TAMANHO_FIGURA = (10, 6)
DPI = 300
PASTA_RESULTADOS = "resultados"
PASTA_FIGURAS = "figuras"
MAX_ITER = 100000
PIVOT_TOL = 1e-12


# -----------------------------------------------------------------------------
# FUNÇÕES AUXILIARES
# -----------------------------------------------------------------------------
def criar_diretorios():
    os.makedirs(PASTA_RESULTADOS, exist_ok=True)
    os.makedirs(PASTA_FIGURAS, exist_ok=True)


def ler_float(texto, minimo=None, maximo=None):
    while True:
        try:
            valor = float(input(texto))
            if minimo is not None and valor <= minimo:
                print(f"Digite um valor maior que {minimo}.")
                continue
            if maximo is not None and valor >= maximo:
                print(f"Digite um valor menor que {maximo}.")
                continue
            return valor
        except ValueError:
            print("Entrada inválida.")


def ler_int(texto, minimo=2):
    while True:
        try:
            valor = int(input(texto))
            if valor < minimo:
                print(f"Digite um inteiro maior ou igual a {minimo}.")
                continue
            return valor
        except ValueError:
            print("Entrada inválida.")


def ler_omega():
    while True:
        omega = ler_float("Fator de relaxação (1 < ω < 2): ")
        if 1.0 < omega < 2.0:
            return omega
        print("Valor inválido.")


def obter_parametros():
    print("\n" + "=" * 70)
    print("CONDUÇÃO DE CALOR EM ALETA RETANGULAR")
    print("=" * 70)
    L = ler_float("Comprimento da aleta (m): ", minimo=0)
    esp = ler_float("Espessura da aleta (m): ", minimo=0)
    k = ler_float("Condutividade térmica (W/m.K): ", minimo=0)
    h = ler_float("Coeficiente convectivo (W/m².K): ", minimo=0)
    Tb = ler_float("Temperatura da base (°C): ")
    Tinf = ler_float("Temperatura ambiente (°C): ")
    nx = ler_int("Número de nós em X: ")
    ny = ler_int("Número de nós em Y: ")
    tolerancia = ler_float("Tolerância: ", minimo=0)
    omega = ler_omega()
    return {
        "L": L, "esp": esp, "k": k, "h": h,
        "Tb": Tb, "Tinf": Tinf,
        "nx": nx, "ny": ny,
        "tol": tolerancia, "omega": omega
    }


def imprimir_parametros(dados):
    print("\n" + "=" * 70)
    print("DADOS DA SIMULAÇÃO")
    print("=" * 70)
    for chave, valor in dados.items():
        print(f"{chave.capitalize():<20}: {valor}")
    print("=" * 70)


def verificar_parametros(dados):
    if dados["Tb"] <= dados["Tinf"]:
        print("\nAviso: Temperatura da base menor ou igual à ambiente.")
    if dados["nx"] < 3 or dados["ny"] < 3:
        raise ValueError("Utilize pelo menos 3 nós em cada direção.")


def inicializar():
    criar_diretorios()
    dados = obter_parametros()
    verificar_parametros(dados)
    imprimir_parametros(dados)
    return dados


# -----------------------------------------------------------------------------
# MALHA E SISTEMA LINEAR
# -----------------------------------------------------------------------------
def construir_malha(dados):
    x = np.linspace(0.0, dados["L"], dados["nx"])
    y = np.linspace(0.0, dados["esp"], dados["ny"])
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X, Y = np.meshgrid(x, y, indexing="ij")
    return x, y, X, Y, dx, dy


def criar_indices(nx, ny):
    indice = np.zeros((nx, ny), dtype=int)
    cont = 0
    for i in range(nx):
        for j in range(ny):
            indice[i, j] = cont
            cont += 1
    return indice


def montar_sistema(dados, dx, dy):
    nx, ny = dados["nx"], dados["ny"]
    k, h = dados["k"], dados["h"]
    Tb, Tinf = dados["Tb"], dados["Tinf"]
    beta_x = h * dx / k
    beta_y = h * dy / k
    N = nx * ny
    A = np.zeros((N, N))
    b = np.zeros(N)
    indice = criar_indices(nx, ny)

    for i in range(nx):
        for j in range(ny):
            p = indice[i, j]

            # Base (Dirichlet)
            if i == 0:
                A[p, p] = 1.0
                b[p] = Tb
                continue

            # Nó interno
            if 0 < i < nx - 1 and 0 < j < ny - 1:
                A[p, p] = 4.0
                A[p, indice[i+1, j]] = -1.0
                A[p, indice[i-1, j]] = -1.0
                A[p, indice[i, j+1]] = -1.0
                A[p, indice[i, j-1]] = -1.0
                continue

            # Extremidade direita (Robin)
            if i == nx - 1 and 0 < j < ny - 1:
                A[p, p] = 4.0 + 2.0 * beta_x
                A[p, indice[i-1, j]] = -2.0
                A[p, indice[i, j+1]] = -1.0
                A[p, indice[i, j-1]] = -1.0
                b[p] = 2.0 * beta_x * Tinf
                continue

            # Borda superior (Robin)
            if j == ny - 1 and 0 < i < nx - 1:
                A[p, p] = 4.0 + 2.0 * beta_y
                A[p, indice[i+1, j]] = -1.0
                A[p, indice[i-1, j]] = -1.0
                A[p, indice[i, j-1]] = -2.0
                b[p] = 2.0 * beta_y * Tinf
                continue

            # Borda inferior (Robin)
            if j == 0 and 0 < i < nx - 1:
                A[p, p] = 4.0 + 2.0 * beta_y
                A[p, indice[i+1, j]] = -1.0
                A[p, indice[i-1, j]] = -1.0
                A[p, indice[i, j+1]] = -2.0
                b[p] = 2.0 * beta_y * Tinf
                continue

            # Canto superior direito
            if i == nx - 1 and j == ny - 1:
                A[p, p] = 4.0 + 2.0 * beta_x + 2.0 * beta_y
                A[p, indice[i-1, j]] = -2.0
                A[p, indice[i, j-1]] = -2.0
                b[p] = (2.0 * beta_x + 2.0 * beta_y) * Tinf
                continue

            # Canto inferior direito
            if i == nx - 1 and j == 0:
                A[p, p] = 4.0 + 2.0 * beta_x + 2.0 * beta_y
                A[p, indice[i-1, j]] = -2.0
                A[p, indice[i, j+1]] = -2.0
                b[p] = (2.0 * beta_x + 2.0 * beta_y) * Tinf
                continue

    return A, b, indice, beta_x, beta_y


def preparar_problema(dados):
    x, y, X, Y, dx, dy = construir_malha(dados)
    A, b, indice, beta_x, beta_y = montar_sistema(dados, dx, dy)
    return x, y, X, Y, dx, dy, A, b, indice, beta_x, beta_y


# -----------------------------------------------------------------------------
# ELIMINAÇÃO DE GAUSS
# -----------------------------------------------------------------------------
def pivotamento_parcial(A, b, k):
    linha = np.argmax(np.abs(A[k:, k])) + k
    if abs(A[linha, k]) < PIVOT_TOL:
        raise ValueError("Matriz singular.")
    if linha != k:
        A[[k, linha]] = A[[linha, k]]
        b[[k, linha]] = b[[linha, k]]


def substituicao_retroativa(A, b):
    n = len(b)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if abs(A[i, i]) < PIVOT_TOL:
            raise ValueError("Divisão por zero.")
        x[i] = (b[i] - np.dot(A[i, i+1:], x[i+1:])) / A[i, i]
    return x


def eliminacao_gauss(A, b):
    A = A.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)
    for k in range(n - 1):
        pivotamento_parcial(A, b, k)
        for i in range(k + 1, n):
            fator = A[i, k] / A[k, k]
            A[i, k:] -= fator * A[k, k:]
            b[i] -= fator * b[k]
    return substituicao_retroativa(A, b)


def calcular_residuo(A, x, b):
    return np.linalg.norm(A @ x - b)


# -----------------------------------------------------------------------------
# SOLUÇÃO ANALÍTICA 1D
# -----------------------------------------------------------------------------
def solucao_analitica(dados, x):
    L, esp, k, h, Tb, Tinf = dados["L"], dados["esp"], dados["k"], dados["h"], dados["Tb"], dados["Tinf"]
    largura = 1.0
    Ac = esp * largura
    P = 2 * (largura + esp)
    m = np.sqrt(h * P / (k * Ac))
    theta = np.zeros_like(x)
    den = np.cosh(m * L) + (h / (m * k)) * np.sinh(m * L)
    for i, xi in enumerate(x):
        num = np.cosh(m * (L - xi)) + (h / (m * k)) * np.sinh(m * (L - xi))
        theta[i] = num / den
    return Tinf + (Tb - Tinf) * theta


def calcular_erros(numerica, analitica):
    erro = np.abs(numerica - analitica)
    return np.mean(erro), np.max(erro), np.sqrt(np.mean(erro**2))


def resolver_gauss(A, b, dados, x, ny):
    inicio = time.time()
    T = eliminacao_gauss(A, b)
    tempo = time.time() - inicio
    residuo = calcular_residuo(A, T, b)
    T2D = T.reshape((dados["nx"], ny))
    perfil = T2D[:, ny // 2]
    analitica = solucao_analitica(dados, x)
    erro_medio, erro_maximo, rmse = calcular_erros(perfil, analitica)

    print("\n" + "=" * 70)
    print("ELIMINAÇÃO DE GAUSS")
    print("=" * 70)
    print(f"Tempo     : {tempo:.6f} s")
    print(f"Resíduo   : {residuo:.6e}")
    print(f"Erro Médio: {erro_medio:.6e}")
    print(f"Erro Máx  : {erro_maximo:.6e}")
    print(f"RMSE      : {rmse:.6e}")
    print("=" * 70)

    return T, analitica, tempo, residuo, erro_medio, erro_maximo, rmse


# -----------------------------------------------------------------------------
# MÉTODO DE LIEBMANN (GAUSS-SEIDEL)
# -----------------------------------------------------------------------------
def criar_chute_inicial(dados, indice):
    T = np.full(dados["nx"] * dados["ny"], dados["Tinf"], dtype=float)
    for j in range(dados["ny"]):
        T[indice[0, j]] = dados["Tb"]
    return T


def atualizar_no(T, indice, i, j, dados, beta_x, beta_y):
    nx, ny = dados["nx"], dados["ny"]
    Tinf = dados["Tinf"]

    if 0 < i < nx-1 and 0 < j < ny-1:
        return (T[indice[i+1, j]] + T[indice[i-1, j]] +
                T[indice[i, j+1]] + T[indice[i, j-1]]) / 4.0

    if i == nx-1 and 0 < j < ny-1:
        return (2*T[indice[i-1, j]] + T[indice[i, j+1]] +
                T[indice[i, j-1]] + 2*beta_x*Tinf) / (4 + 2*beta_x)

    if j == ny-1 and 0 < i < nx-1:
        return (T[indice[i+1, j]] + T[indice[i-1, j]] +
                2*T[indice[i, j-1]] + 2*beta_y*Tinf) / (4 + 2*beta_y)

    if j == 0 and 0 < i < nx-1:
        return (T[indice[i+1, j]] + T[indice[i-1, j]] +
                2*T[indice[i, j+1]] + 2*beta_y*Tinf) / (4 + 2*beta_y)

    if i == nx-1 and j == ny-1:
        return (2*T[indice[i-1, j]] + 2*T[indice[i, j-1]] +
                2*(beta_x + beta_y)*Tinf) / (4 + 2*(beta_x + beta_y))

    if i == nx-1 and j == 0:
        return (2*T[indice[i-1, j]] + 2*T[indice[i, j+1]] +
                2*(beta_x + beta_y)*Tinf) / (4 + 2*(beta_x + beta_y))

    return T[indice[i, j]]


def metodo_liebmann(dados, indice, beta_x, beta_y):
    T = criar_chute_inicial(dados, indice)
    historico = []
    erro = np.inf
    iteracoes = 0
    inicio = time.time()

    while erro > dados["tol"] and iteracoes < MAX_ITER:
        erro = 0.0
        for i in range(1, dados["nx"]):
            for j in range(dados["ny"]):
                p = indice[i, j]
                antigo = T[p]
                novo = atualizar_no(T, indice, i, j, dados, beta_x, beta_y)
                T[p] = novo
                erro = max(erro, abs(novo - antigo))
        historico.append(erro)
        iteracoes += 1

    tempo = time.time() - inicio
    convergiu = erro <= dados["tol"]
    return T, historico, iteracoes, erro, tempo, convergiu


def imprimir_liebmann(tempo, iteracoes, erro, convergiu):
    print("\n" + "=" * 70)
    print("MÉTODO DE LIEBMANN")
    print("=" * 70)
    print(f"Tempo       : {tempo:.6f} s")
    print(f"Iterações   : {iteracoes}")
    print(f"Erro Final  : {erro:.6e}")
    print(f"Convergência: {convergiu}")
    print("=" * 70)


# -----------------------------------------------------------------------------
# MÉTODO SOR
# -----------------------------------------------------------------------------
def metodo_sor(dados, indice, beta_x, beta_y):
    T = criar_chute_inicial(dados, indice)
    historico = []
    erro = np.inf
    iteracoes = 0
    omega = dados["omega"]
    inicio = time.time()

    while erro > dados["tol"] and iteracoes < MAX_ITER:
        erro = 0.0
        for i in range(1, dados["nx"]):
            for j in range(dados["ny"]):
                p = indice[i, j]
                antigo = T[p]
                valor_gs = atualizar_no(T, indice, i, j, dados, beta_x, beta_y)
                novo = antigo + omega * (valor_gs - antigo)
                T[p] = novo
                erro = max(erro, abs(novo - antigo))
        historico.append(erro)
        iteracoes += 1

    tempo = time.time() - inicio
    convergiu = erro <= dados["tol"]
    return T, historico, iteracoes, erro, tempo, convergiu


def imprimir_sor(tempo, iteracoes, erro, convergiu):
    print("\n" + "=" * 70)
    print("MÉTODO SOR")
    print("=" * 70)
    print(f"Tempo       : {tempo:.6f} s")
    print(f"Iterações   : {iteracoes}")
    print(f"Erro Final  : {erro:.6e}")
    print(f"Convergência: {convergiu}")
    print("=" * 70)


# -----------------------------------------------------------------------------
# COMPARAÇÕES E ESTUDOS
# -----------------------------------------------------------------------------
def comparar_metodos(dados, T_gauss, T_liebmann, T_sor,
                     tempo_gauss, tempo_liebmann, tempo_sor,
                     iter_liebmann, iter_sor):
    print("\n" + "=" * 80)
    print("COMPARAÇÃO DOS MÉTODOS")
    print("=" * 80)
    print(f"{'Método':<20}{'Tempo(s)':>15}{'Iterações':>18}")
    print("-" * 80)
    print(f"{'Gauss':<20}{tempo_gauss:>15.6f}{'-':>18}")
    print(f"{'Liebmann':<20}{tempo_liebmann:>15.6f}{iter_liebmann:>18}")
    print(f"{'SOR':<20}{tempo_sor:>15.6f}{iter_sor:>18}")
    print("=" * 80)

    erro_GL = np.max(np.abs(T_gauss - T_liebmann))
    erro_GS = np.max(np.abs(T_gauss - T_sor))
    erro_LS = np.max(np.abs(T_liebmann - T_sor))
    print("\nDiferença Máxima entre Soluções")
    print(f"Gauss  x Liebmann : {erro_GL:.6e}")
    print(f"Gauss  x SOR      : {erro_GS:.6e}")
    print(f"Liebmann x SOR    : {erro_LS:.6e}")
    print("=" * 80)


def estudo_omega(dados, indice, beta_x, beta_y):
    print("\n" + "=" * 70)
    print("ESTUDO DO FATOR DE RELAXAÇÃO")
    print("=" * 70)
    omegas = np.arange(1.0, 2.0, 0.1)
    iteracoes = []
    for omega in omegas:
        dados_temp = dados.copy()
        dados_temp["omega"] = omega
        _, _, it, erro, _, _ = metodo_sor(dados_temp, indice, beta_x, beta_y)
        iteracoes.append(it)
        print(f"ω = {omega:.1f}   Iterações = {it:5d}   Erro = {erro:.2e}")
    return omegas, iteracoes


def melhor_omega(omegas, iteracoes):
    idx = np.argmin(iteracoes)
    print("\n" + "=" * 70)
    print("MELHOR FATOR ω")
    print("=" * 70)
    print(f"ω = {omegas[idx]:.1f}  (Iterações = {iteracoes[idx]})")
    print("=" * 70)


# -----------------------------------------------------------------------------
# GRÁFICOS E EXPORTAÇÃO
# -----------------------------------------------------------------------------
def salvar_resultados(dados, x, y, indice, T_gauss, T_liebmann, T_sor, T_analitica):
    arquivo = os.path.join(PASTA_RESULTADOS, "dados_simulacao.csv")
    with open(arquivo, "w", encoding="utf-8") as f:
        f.write("x;y;Gauss;Liebmann;SOR;Analitica\n")
        for i in range(dados["nx"]):
            for j in range(dados["ny"]):
                p = indice[i, j]
                f.write(f"{x[i]:.6f};{y[j]:.6f};"
                        f"{T_gauss[p]:.6f};{T_liebmann[p]:.6f};"
                        f"{T_sor[p]:.6f};{T_analitica[i]:.6f}\n")
    print(f"\nArquivo salvo em {arquivo}")


def grafico_temperatura(X, Y, T2D):
    plt.figure(figsize=TAMANHO_FIGURA)
    plt.contourf(X, Y, T2D, 20, cmap="hot")
    plt.colorbar(label="Temperatura (°C)")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Mapa de Temperatura")
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_FIGURAS, "mapa_temperatura.png"), dpi=DPI)
    plt.show()


def grafico_perfil(dados, x, T_gauss, T_liebmann, T_sor, T_analitica):
    ny = dados["ny"]
    linha = ny // 2
    plt.figure(figsize=TAMANHO_FIGURA)
    plt.plot(x, T_gauss.reshape((-1, ny))[:, linha], label="Gauss")
    plt.plot(x, T_liebmann.reshape((-1, ny))[:, linha], "--", label="Liebmann")
    plt.plot(x, T_sor.reshape((-1, ny))[:, linha], "-.", label="SOR")
    plt.plot(x, T_analitica, ":", linewidth=3, label="Analítica")
    plt.xlabel("x (m)")
    plt.ylabel("Temperatura (°C)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_FIGURAS, "perfil_central.png"), dpi=DPI)
    plt.show()


def grafico_convergencia(hist_liebmann, hist_sor):
    plt.figure(figsize=TAMANHO_FIGURA)
    plt.semilogy(hist_liebmann, label="Liebmann")
    plt.semilogy(hist_sor, label="SOR")
    plt.xlabel("Iterações")
    plt.ylabel("Erro")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_FIGURAS, "convergencia.png"), dpi=DPI)
    plt.show()


def estudo_refinamento(dados):
    print("\n" + "=" * 70)
    print("ESTUDO DE REFINAMENTO")
    print("=" * 70)
    malhas = [(10, 5), (20, 10), (30, 15), (40, 20), (60, 30)]
    erros = []
    elementos = []
    for nx, ny in malhas:
        dados_temp = dados.copy()
        dados_temp["nx"] = nx
        dados_temp["ny"] = ny
        x, _, _, _, _, _, A, b, _, _, _ = preparar_problema(dados_temp)
        T = eliminacao_gauss(A, b)
        T2D = T.reshape((nx, ny))
        perfil = T2D[:, ny // 2]
        analitica = solucao_analitica(dados_temp, x)
        erro = np.mean(np.abs(perfil - analitica))
        elementos.append(nx * ny)
        erros.append(erro)
        print(f"{nx:3d} x {ny:<3d} -> erro médio = {erro:.6e}")

    plt.figure(figsize=TAMANHO_FIGURA)
    plt.loglog(elementos, erros, "o-")
    plt.xlabel("Número de Nós")
    plt.ylabel("Erro Médio")
    plt.title("Refinamento da Malha")
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_FIGURAS, "refinamento.png"), dpi=DPI)
    plt.show()
    return elementos, erros


# -----------------------------------------------------------------------------
# FUNÇÃO PRINCIPAL
# -----------------------------------------------------------------------------
def main():
    print("\n" + "=" * 70)
    print("PPC6 - MÉTODO DAS DIFERENÇAS FINITAS")
    print("=" * 70)

    dados = inicializar()
    x, y, X, Y, dx, dy, A, b, indice, beta_x, beta_y = preparar_problema(dados)

    # Gauss
    T_gauss, T_analitica, tempo_gauss, residuo, erro_medio, erro_maximo, rmse = \
        resolver_gauss(A, b, dados, x, dados["ny"])

    # Liebmann
    T_liebmann, hist_liebmann, iter_liebmann, erro_liebmann, tempo_liebmann, conv_lieb = \
        metodo_liebmann(dados, indice, beta_x, beta_y)
    imprimir_liebmann(tempo_liebmann, iter_liebmann, erro_liebmann, conv_lieb)

    # SOR
    T_sor, hist_sor, iter_sor, erro_sor, tempo_sor, conv_sor = \
        metodo_sor(dados, indice, beta_x, beta_y)
    imprimir_sor(tempo_sor, iter_sor, erro_sor, conv_sor)

    # Comparações
    comparar_metodos(dados, T_gauss, T_liebmann, T_sor,
                     tempo_gauss, tempo_liebmann, tempo_sor,
                     iter_liebmann, iter_sor)

    # Estudo ω
    omegas, iteracoes = estudo_omega(dados, indice, beta_x, beta_y)
    melhor_omega(omegas, iteracoes)

    # Gráficos
    grafico_temperatura(X, Y, T_gauss.reshape((dados["nx"], dados["ny"])))
    grafico_perfil(dados, x, T_gauss, T_liebmann, T_sor, T_analitica)
    grafico_convergencia(hist_liebmann, hist_sor)

    # Exportação
    salvar_resultados(dados, x, y, indice, T_gauss, T_liebmann, T_sor, T_analitica)

    # Refinamento
    estudo_refinamento(dados)

    # Resumo final
    print("\n" + "=" * 70)
    print("RESUMO FINAL")
    print("=" * 70)
    print(f"Nós da malha      : {dados['nx']} x {dados['ny']}")
    print(f"Tempo Gauss       : {tempo_gauss:.6f} s")
    print(f"Tempo Liebmann    : {tempo_liebmann:.6f} s")
    print(f"Tempo SOR         : {tempo_sor:.6f} s")
    print(f"Erro Médio        : {erro_medio:.6e}")
    print(f"Erro Máximo       : {erro_maximo:.6e}")
    print(f"RMSE              : {rmse:.6e}")
    print(f"Resíduo           : {residuo:.6e}")
    print("=" * 70)
    print("\nSimulação concluída com sucesso!")
    print("Arquivos gerados:")
    print("   resultados/dados_simulacao.csv")
    print("   figuras/mapa_temperatura.png")
    print("   figuras/perfil_central.png")
    print("   figuras/convergencia.png")
    print("   figuras/refinamento.png")


# -----------------------------------------------------------------------------
# EXECUÇÃO
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
