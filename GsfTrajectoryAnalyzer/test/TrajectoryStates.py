class SingleTState:
    def __init__(self,parameters,errors,weight=1.):
        self.pars_ = parameters
        self.errors_ = errors
        self.weight_ = weight

    def weight(self):
        return self.weight_

    def __getitem__(self,ind):
        return (self.pars_[ind], self.errs_[ind])

class MultiTState:

    def __init__(self,parameters=None,errors=None):
        self.meanPars_ = parameters
        self.meanErrs_ = errors
        self.nc_ = 0
        self.components_ = [ ]

    def addComponent(self,parameters,errors,weight):
        self.nc_ += 1
        self.components_.append(SingleTState(parameters,errors,weight))

    def weights(self):
        return [ c.weight() for c in self.components_ ]

    def __getitem__(self,ind):
        return ( [ c[ind][0] for c in self.components_ ], [ c[ind][1] for c in self.components_ ] )

