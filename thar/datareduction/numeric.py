def f(x):
    return x**10+x**9+x**8+x**7+x**6+x**5+x**4+x**3+x**2+x+1

def ff(x):
    return f(x) - (x**11-1) / (x-1)