import json
from sage.all import ZZ, valuation, product, factor

class CS_Params:

    def __init__(self, a=None, filename=None, dic=None):
        """
        Helper class for parameters. Can be loaded from
        a dictionary or a json file.
        """
        if a:
            filename = f'params/cs_ftf{a}.json'
        if filename:
            self.load_from_file(filename)
        elif dic:
            self.load_from_dic(dic)
        else:
            raise ValueError('provide either a filename or a dictionary')

    def load_from_dic(self, dic):
        self.a = ZZ(dic['a'])
        self.B = ZZ(dic['B'])
        self.p = ZZ(dic['p'])
        self.cof = ZZ(dic['cof']) # p = 2^mu*[primes]*cof + 1

        self.primes = [ZZ(mi) for mi in dic['primes']]

        self.l1 = ZZ(dic['l1'])
        self.l2 = ZZ(dic['l2'])
        self.alpha = ZZ(dic['alpha'])
        self.beta = ZZ(dic['beta'])
        self.eigenvals = {ZZ(i):j for i,j in dic['eigenvals'].items()}

        # Order of the 2 torsion and of the rational torsion to use
        self.eval_even = self.a+1
        self.eval_odd = product([i for (i, _) in factor(self.p+1)[1:]])

    def load_from_file(self, filename):
        dic = json.load(open(filename, 'r'))
        self.load_from_dic(dic)

class sub_Params:

    def __init__(self, a=None, filename=None, dic=None):
        """
        Helper class for parameters. Can be loaded from
        a dictionary or a json file.
        """
        if a:
            filename = f'params/sub_pgen{a}.json'
        if filename:
            self.load_from_file(filename)
        elif dic:
            self.load_from_dic(dic)
        else:
            raise ValueError('provide either a filename or a dictionary')

    def load_from_dic(self, dic):
        self.a = ZZ(dic['a'])
        self.B = ZZ(dic['B'])
        self.p = ZZ(dic['p'])
        self.cof = ZZ(dic['cof'])
        self.primes = [ZZ(mi) for mi in dic['primes']]

        self.f1 = ZZ(dic['f1'])
        self.f2 = ZZ(dic['f2'])
        self.f = ZZ(dic['f'])
        self.eval_even = ZZ(dic['eval_even'])
        self.eval_odd = ZZ(dic['eval_odd'])
        self.m = ZZ(dic['m'])
        self.n = ZZ(dic['n'])

        self.eigenvals = {ZZ(i):j for i,j in dic['eigenvals'].items()}

        self.E = [[ZZ(dic['E'][i][0]), ZZ(dic['E'][i][1])] for i in range(0,len(dic['E']))] 
        self.E_tors = [ZZ(i) for i in dic['E_tors']]
        self.tau = [ZZ(i) for i in dic['tau']] 

    def load_from_file(self, filename):
        dic = json.load(open(filename, 'r'))
        self.load_from_dic(dic)


