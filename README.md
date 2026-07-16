# PPC6 - Método das Diferenças Finitas para Condução de Calor em Aletas

## Disciplina

Cálculo Numérico Aplicado

Universidade de Brasília (UnB)

---

## Descrição

Este projeto implementa a solução numérica da distribuição de temperatura em uma aleta retangular utilizando o Método das Diferenças Finitas (MDF).

O programa compara três métodos de solução:

- Eliminação de Gauss
- Método de Liebmann (Gauss-Seidel)
- Método SOR (Successive Over-Relaxation)

Além disso, o programa compara os resultados numéricos com a solução analítica, realiza um estudo de refinamento da malha, gera gráficos e exporta os resultados para arquivo CSV.

---

# Funcionalidades

- Construção automática da malha
- Montagem do sistema linear
- Eliminação de Gauss com pivotamento parcial
- Método de Liebmann
- Método SOR
- Solução analítica da aleta
- Comparação entre métodos
- Estudo do fator de relaxação (ω)
- Estudo de refinamento da malha
- Geração automática de gráficos
- Exportação dos resultados em CSV

---

# Estrutura do projeto

```
PPC6/
│
├── PPC6.py
├── README.md
├── resultados/
│   └── dados_simulacao.csv
│
└── figuras/
    ├── mapa_temperatura.png
    ├── perfil_central.png
    ├── convergencia.png
    └── refinamento.png
```

---

# Bibliotecas utilizadas

O projeto utiliza apenas bibliotecas padrão do Python e NumPy/Matplotlib.

```python
numpy
matplotlib
time
os
```

---

# Como executar

## 1. Instale o Python

Recomendado:

Python 3.10 ou superior

---

## 2. Instale as bibliotecas

Abra o terminal e execute:

```bash
pip install numpy matplotlib
```

ou

```bash
python -m pip install numpy matplotlib
```

---

## 3. Execute o programa

No terminal:

```bash
python PPC6.py
```

Caso tenha mais de uma versão do Python:

```bash
python3 PPC6.py
```

---

# Dados solicitados pelo programa

Durante a execução serão solicitados:

- Comprimento da aleta
- Espessura
- Condutividade térmica
- Coeficiente convectivo
- Temperatura da base
- Temperatura ambiente
- Número de nós em X
- Número de nós em Y
- Tolerância
- Fator de relaxação (ω)

---

# Métodos implementados

## Eliminação de Gauss

Método direto para solução do sistema linear.

Características:

- Pivotamento parcial
- Alta precisão
- Utilizado como referência

---

## Método de Liebmann

Método iterativo baseado em Gauss-Seidel.

Características:

- Atualização sucessiva dos nós
- Controle de tolerância
- Histórico de convergência

---

## Método SOR

Método iterativo utilizando sobre-relaxação.

Características:

- Utiliza o fator ω
- Geralmente converge mais rápido que Gauss-Seidel
- Permite estudo do melhor fator de relaxação

---

# Saídas do programa

Ao final da execução são gerados:

## Arquivos

```
resultados/dados_simulacao.csv
```

## Figuras

```
mapa_temperatura.png

perfil_central.png

convergencia.png

refinamento.png
```

---

# Resultados apresentados

O programa exibe:

- Tempo de execução
- Número de iterações
- Resíduo
- Erro médio
- Erro máximo
- RMSE
- Comparação entre os métodos
- Melhor fator ω
- Refinamento da malha

---

# Organização do código

O código foi dividido em módulos lógicos:

1. Entrada de dados
2. Construção da malha
3. Montagem do sistema
4. Eliminação de Gauss
5. Solução analítica
6. Método de Liebmann
7. Método SOR
8. Comparação entre métodos
9. Gráficos
10. Refinamento da malha
11. Função principal

---

# Observações

O projeto foi desenvolvido para fins acadêmicos na disciplina de Cálculo Numérico Aplicado.

