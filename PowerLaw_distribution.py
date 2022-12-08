# 为了保证x的取值范围是[epsilon,1],确保概率密度函数的积分为1
def normalizing_constant(gamma, epsilon):
    return (1 - epsilon**(1 - gamma)) / (1 - gamma)

# 概率密度函数
def pdf(gamma, epsilon):
    normalizing_const = normalizing_constant(gamma, epsilon)
    def normalized_pdf(x):
        if not epsilon <= x <= 1.0:
            return 0
        return x**(-gamma) / normalizing_const
    return normalized_pdf

# 累积分布函数
def cdf(gamma, epsilon):
    normalizing_const = normalizing_constant(gamma, epsilon)
    def normalized_cdf(x):
        if x < epsilon:
            return 0
        elif x >= 1:
            return 1
        a = x**(1 - gamma) - epsilon**(1 - gamma)
        b = normalizing_const * (1 - gamma)
        return a / b
    return normalized_cdf

# 逆变换抽样
def inv_cdf(gamma, epsilon):
    normalizing_const = normalizing_constant(gamma, epsilon)
    def normalized_inv_cdf(x):
        if not 0 <= x <= 1:
            raise ValueError()
        a = x * normalizing_const * (1 - gamma)
        b = epsilon**(1 - gamma)
        return (a + b)**(1 / (1 - gamma))
    return normalized_inv_cdf

# 求期望/一阶原点矩
# def expected_value(gamma, epsilon):
#     normalizing_const = normalizing_constant(gamma, epsilon)
#     return (1.0 - epsilon**(2.0 - gamma)) / (normalizing_const * (2 - gamma))

# E = expected_value(2.8, 0.01)
# print(E)

#求二阶原点矩
# def second_moment(gamma, epsilon):
#     normalizing_const = normalizing_constant(gamma, epsilon)
#     return (1.0 - epsilon**(3.0 - gamma))/(normalizing_const * (3 - gamma))

# E2 = second_moment(2.8, 0.01)
# print(E2)

# 抽样
# def inverse_transform_sampling(inv_cdf):
#     x = inv_cdf
#     while True:
#         yield x


#带yield的函数是一个生成器，而不是一个函数了，这个生成器有一个函数就是next函数，
# next就相当于“下一步”生成哪个数，这一次的next开始的地方是接着上一次的next停止的地方执行的，
# 所以调用next的时候，生成器并不会从foo函数的开始执行，只是接着上一步停止的地方开始，
# 然后遇到yield后，return出要生成的数，此步就结束。




#### 固定上下限
# 为了保证m的取值范围是[m_min,m_max],确保概率密度函数的积分为1
def m_normalizing_constant(gamma, m_min, m_max):
    return (m_max**(1-gamma) - m_min**(1 - gamma)) / (1 - gamma)

# 概率密度函数
def m_pdf(gamma, m_min, m_max):
    m_normalizing_const = m_normalizing_constant(gamma, m_min, m_max)
    def m_normalized_pdf(x):
        if not m_min <= x <= m_max:
            return 0
        return x**(-gamma) / m_normalizing_const
    return m_normalized_pdf

# 累积分布函数
def m_cdf(gamma, m_min, m_max):
    m_normalizing_const = m_normalizing_constant(gamma, m_min, m_max)
    def m_normalized_cdf(x):
        if x < m_min:
            return 0
        elif x >= m_max:
            return 1
        a = x**(1 - gamma) - m_min**(1 - gamma)
        b = m_normalizing_const * (1 - gamma)
        return a / b
    return m_normalized_cdf

# 逆变换抽样
def m_inv_cdf(gamma, m_min, m_max):
    m_normalizing_const = m_normalizing_constant(gamma, m_min, m_max)
    def m_normalized_inv_cdf(x):
        if not 0 <= x <= 1:
            raise ValueError()
        a = x * m_normalizing_const * (1 - gamma)
        b = m_min**(1 - gamma)
        return (a + b)**(1 / (1 - gamma))
    return m_normalized_inv_cdf


# # 求期望/一阶原点矩
def expected_value(gamma, m_min, m_max):
    m_normalizing_const = m_normalizing_constant(gamma, m_min, m_max)
    return (m_max**(2.0-gamma) - m_min**(2.0 - gamma)) / (m_normalizing_const * (2 - gamma))

E = expected_value(2.8, 1, 10)
print(E)

# #求二阶原点矩
# def second_moment(gamma, m_min, m_max):
#     m_normalizing_const = m_normalizing_constant(gamma, m_min, m_max)
#     return (m_max**(3.0-gamma) - m_min**(3.0 - gamma))/(m_normalizing_const * (3 - gamma))
#
# E2 = second_moment(2.8, 1, 10)
# print(E2)


























