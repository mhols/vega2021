import numpy as np
__doc__="linear mixed models"

class LinearMixedModel:
    """
    general linear mixed model
    """

    def __init__(self, d, F,  S=1.0, Z=None, G=None, iG=None):
        self.d = d
        self.F = F
        self.S = S
        self.Z = Z
        self.G = G
        self.iG = iG
        self.LS = None # Cholesky factor of S

        if self.nr == 0:
            self._posterior_lm()
        else:
            self._posterior_lmm()

    def _prewhitening(self, A):

        if type(self.S) is float:
            return A / np.sqrt(self.S)
        if len(self.S) ==  self.F.shape[0]:
            return np.sqrt(1./self.S)[:, np.newaxis] * A
        if self.LS is None:
            self. LS = np.linalg.cholesky(self.S)
            return np.linalg.solve(self.LS, A)
        return A

    @property
    def nr(self):
        if self.Z is None:
            return 0
        else:
            return self.Z.shape[1]

    @property
    def nf(self):
        if self.F is None:
            return 0
        else:
            return self.F.shape[1]

    def _posterior_lm(self):
        d = self._prewhitening(self.d)
        F = self._prewhitening(self.F)

        self._posterior_mean = np.linalg.lstsq(F, d)[0]
        self._posterior_inv_variance = np.dot(F.T, F)

    def _posterior_lmm(self):
        d = self._prewhitening(self.d)
        F = self._prewhitening(self.F)
        Z = self._prewhitening(self.Z)

        W = np.column_stack((F,Z))
        X = np.dot(W.T, W)
        if self.iG is None:
            X[-self.nr:, -self.nr:] += np.linalg.inv(self.G)
        else:
            X[-self.nr:, -self.nr:] += self.iG
        self._posterior_inv_variance =  X
        self._posterior_mean = np.linalg.solve(X, np.dot(W.T, d))

    def post_mean_fixed_effects(self):
        return self._posterior_mean[:self.nf]

    def post_mean_random_effects(self):
        return self._posterior_mean[self.nf:]

if __name__ == '__main__':
    
    def K(x,y):
        return np.abs(x-y)**3/12

    x = np.array([1, 1.1, 1.3, 1.6, 2, 3.5, 4.7, 5])

    d = x * (5 - x) + 2 * np.random.normal(size=x.shape[0])

    Z = K(x[:, np.newaxis], x[np.newaxis, :])
    print(Z)
    F = np.zeros((d.shape[0], 3))
    F[:, 0] = 1
    F[:, 1] = x
    F[:, 2] = x**2

    myModel = LinearMixedModel(d, F, 1.0, Z=np.identity(x.shape[0]), G = Z)
    print(myModel._posterior_mean)

