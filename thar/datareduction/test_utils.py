import unittest
import util
import numpy as np

class TestUtil(unittest.TestCase):

    def test_translation_same_grid(self):

        def f(x):
            return np.sin((x/50)**2)/ (1+(x/100)**2)

        d = 24.78
        def ff(x):
            return 3.5*f(x - d)
        
        N = 1000

        x = np.arange(N)
        y = f(x) + 0.1 * np.random.uniform(size=x.shape[0])


        xx = 3.2 + np.arange(N/2)

        yy = ff(xx) +  0.1 * np.random.uniform(size=xx.shape[0])

        dd = util.translation_same_grid(x[0],y, xx[0], yy, d-10,d+10)
        

        self.assertTrue( np.abs(d-dd) < 0.7)
    
    def test_translation_same_grid_2(self):

        N = 400
        M = 55
        a = 48
        aa = 23

        p =44
        pp = -2334

        x = np.zeros(N)
        xx = np.zeros(M)

        x[a] = 1
        x[a+1] =1
        xx[aa] = 1
        xx[aa+1] = 1

        d = aa + pp - (a + p)
        dd = util.translation_same_grid(p, x, pp, xx, d-20)
        self.assertTrue( d == dd )
        dd = util.translation_same_grid(p, x, pp, xx, None, d+20)
        self.assertTrue( d == dd )
        dd = util.translation_same_grid(p, x, pp, xx, d-20, d)
        self.assertTrue( d == dd )
        dd = util.translation_same_grid(p, x, pp, xx, d, d)
        self.assertTrue( d == dd )
        dd = util.translation_same_grid(p, x, pp, xx, d-20, d+2000)
        self.assertTrue( d == dd )

        #dd = util.translation_same_grid(p, x, pp, xx, d-20, d+2000, full=True)
        #self.assertTrue( np.abs(d -dd) <= 0.1 )

"""
    def test_translation(self):
        d = 0.7

        def f(x):
            return np.sin(x**2)/ (1+x**2)

        def ff(x):
            return f(x-d)

        x = np.linspace(-1.3, 7.2, 1000)
        y = f(x) + 0.001 * np.random.uniform(size=x.shape)

        xx = np.linspace(-1, 9, 2031)

        yy = ff(xx) +  0.001 * np.random.uniform(size=xx.shape)

        dd = util.translation(x,y,xx,yy)
        
        print(d, dd)

        self.assertTrue( np.abs(d-dd) < 0.1)

    
    def test_wasserstein(self):

        a = 2.3
        b = -0.5

        print('truth ',  b, a)
        x = np.linspace(1,4, 1000)
        y = 4 + np.sin(x*x) + 0.1*np.random.uniform(size=1000)

        xx = (x - b) / a
        yy = 2*y + 0.1 * np.random.uniform(size=1000)

        # subset of xx, yy
        xx = xx[60:-100]
        yy = yy[60:-100]

        bb, aa = util.homothetie_wasserstein(x,y,xx,yy, b-0.5, b+0.5, a-0.02, a+0.02)
        print (bb, aa)
        

        self.assertTrue(np.abs(a - aa) < 0.2)
        self.assertTrue(np.abs(b - bb) < 0.2)
"""
if __name__ == '__main__':
    unittest.main()

