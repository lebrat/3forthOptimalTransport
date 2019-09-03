import numpy as np
import ot3D
import matplotlib.pyplot as plt
import scipy.sparse as sps

class LaguerreTri():
    def __init__(self):
        self.__Tri = ot3D.OTproblem()
    def printCGALVERSION(self) :
        self.__Tri.printCGALVERSION()
    def setCoordinates(self,position):
        self.__Tri.setLaguerrePosition(position[:,0],position[:,1],position[:,2])
    def setPolyLine(self,polyLine,weights):
        self.nb_polyLinept=polyLine.shape[0]
        self.__Tri.setPolyline(polyLine[:,0],polyLine[:,1],polyLine[:,2],weights)

    def setPolyLineNew(self,polyLine,weights):
        self.nb_polyLinept=polyLine.shape[0]//2
        self.__Tri.setPolylineNew(polyLine[:,0],polyLine[:,1],polyLine[:,2],weights)

    def getPolyline(self,i):
        
        
        return np.array(self.__Tri.getContributorLine(i))
    def compute(self,psi,Hessian=False):
        self.__Tri.setLaguerrePsi(psi)
        n = psi.shape[0]
        self.n = psi.shape[0]
        mass = np.zeros(n)
        cost = np.zeros(n)
        Bar = np.zeros(3*n)
        self.__Tri.computeAll(mass,cost,Bar)
        Bar = np.reshape(Bar,(n,3),order='C')
        grad = 1./float(n) - mass
        if Hessian :
            R = self.__Tri.computeHessian()
            Rsparse = sps.coo_matrix((np.array(R.val),(np.array(R.row),np.array(R.col))),shape=(R.getDimRow(),R.getDimCol()))
            return mass,Bar,cost,grad,Rsparse
        else :
            return mass,Bar,cost,grad
    def plotEdges(self,ax,i,dim=3) :
        edges=self.__Tri.getadjEdge(i)
        if dim==3 :
            for e in edges :
                ax.plot([e[0],e[3]],[e[1],e[4]],[e[2],e[5]])
        else :
            for e in edges :
                ax.plot([e[0],e[3]],[e[1],e[4]])
    def plotList(self,l,i,dim=3) :
        edges=self.__Tri.getadjEdge(i)
        if dim==3 :
            for e in edges :
                l.append([[e[0],e[3]],[e[1],e[4]],[e[2],e[5]]])
        else :
            for e in edges :
                l.append([[e[0],e[3]],[e[1],e[4]]])
        return l

    def plotIntersection(self,ax,i,dim=3) :
        inter=self.__Tri.getCutPosition(i)
        if dim==3 :
            for e in inter :
                ax.scatter(e[0],e[1],e[2],c='black')
        else :
            for e in inter :
                ax.scatter(e[0],e[1],c='black')
    def derivativePol(self) :
        dmass = np.zeros(6*self.nb_polyLinept)
        dcost = np.zeros(6*self.nb_polyLinept)
        self.__Tri.computeDeriv(dmass,dcost)
        dmass = np.reshape(dmass,(self.nb_polyLinept,6),order='C')
        dcost = np.reshape(dcost,(self.nb_polyLinept,6),order='C')
        return dmass,dcost

    def derivativePolTotal(self) :
        dmass = np.zeros(6*self.nb_polyLinept)
        dcost = np.zeros(6*self.nb_polyLinept)
        drho = np.zeros(6*self.nb_polyLinept)
        self.__Tri.computeDerivTotal(dmass,dcost,drho)
        dmass = np.reshape(dmass,(self.nb_polyLinept,6),order='C')
        dcost = np.reshape(dcost,(self.nb_polyLinept,6),order='C')
        drho  = np.reshape(drho ,(self.nb_polyLinept,6),order='C')
        return dmass,dcost,drho



    def computePolyinfo(self) :
        massSeen = np.zeros(self.nb_polyLinept)
        costPaid = np.zeros(self.nb_polyLinept)
        self.__Tri.computePolyInfo(massSeen,costPaid)

        return massSeen, costPaid

    def setMassDiracs(self,w):
        self.__Tri.setMassDiracs(w)
        
    def perfomOptimisationBFGS(self,gradTol,nMaxIter,w1,w2,wMaxIter,memSize=15,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBFGS(gradTol,nMaxIter,w1,w2,wMaxIter,memSize,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()


    def perfomOptimisationNew(self,gradTol,nMaxIter,w1,w2,wMaxIter,memSize=15,leven=1e-5,psi=None,parallelism=False,nbThreads=-1):        
        if psi is not None:
            self.__Tri.perfomOptimPsiBFGSUnleashDaBeastRestart(psi,gradTol,nMaxIter,w1,w2,wMaxIter,memSize,leven,parallelism,nbThreads)
        else:
            self.__Tri.perfomOptimPsiBFGSUnleashDaBeast(gradTol,nMaxIter,w1,w2,wMaxIter,memSize,leven,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationGradient(self,stepSize,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiConstantGD(stepSize,gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationGradientGreedy(self,steps,ranking,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiGreedyGD(steps,ranking,gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationGradientWolfe(self,stepSize,gradTol,nMaxIter,w1,w2,wMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiStrongWolfGD(stepSize,gradTol,nMaxIter,w1,w2,wMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()
    
    def perfomOptimisationBB(self,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiB(gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationBBG(self,gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiBWolf(gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationBBG2(self,gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,maxIncreaseHidden,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiBWolf2(gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,maxIncreaseHidden,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationNesterov(self,Lips,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBall(Lips,gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationNesterovFred(self,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBallFred(gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationNesterovGreedy(self,steps,ranking,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBallGreedy(steps,ranking,gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()

    def perfomOptimisationPolyakGreedy(self,steps,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiPolyak(steps,gradTol,nMaxIter,parallelism,nbThreads)
        psiOpt = np.zeros(self.n)
        self.__Tri.getOptimizedPsi(psiOpt)
        return psiOpt,self.__Tri.getOptimizationCF(),self.__Tri.getOptimizationCV(),self.__Tri.getOptimizationtime()


#######################################################################################################################
######################   ____  _____ ____  _____    _    ____   ____ _   _   ############################################
######################  |  _ \| ____/ ___|| ____|  / \  |  _ \ / ___| | | |  ############################################
######################  | |_) |  _| \___ \|  _|   / _ \ | |_) | |   | |_| |  ############################################
######################  |  _ <| |___ ___) | |___ / ___ \|  _ <| |___|  _  |  ############################################
######################  |_| \_\_____|____/|_____/_/   \_\_| \_\\____|_| |_|  ############################################
######################                                                       ############################################
#######################################################################################################################

# meme fonction mais avec les historiques moins verbeux

    def perfomOptimisationBFGSResearch(self,gradTol,nMaxIter,w1,w2,wMaxIter,memSize=15,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBFGSResearch(gradTol,nMaxIter,w1,w2,wMaxIter,memSize,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationGradientResearch(self,stepSize,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiConstantGDResearch(stepSize,gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationGradientGreedyResearch(self,steps,ranking,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiGreedyGDResearch(steps,ranking,gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationGradientWolfeResearch(self,stepSize,gradTol,nMaxIter,w1,w2,wMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiStrongWolfGDResearch(stepSize,gradTol,nMaxIter,w1,w2,wMaxIter,parallelism,nbThreads)

        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    
    def perfomOptimisationBBResearch(self,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiBResearch(gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationBBGResearch(self,gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiBWolfResearch(gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime
        
    def perfomOptimisationBBG2Research(self,gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,maxIncreaseHidden,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBarzilaiBWolf2Research(gradTol,lineSearchDecrease,gammaWolfe,nMaxIter,memSize,maxIncreaseHidden,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationNesterovResearch(self,Lips,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBallResearch(Lips,gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationNesterovFredResearch(self,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBallFredResearch(gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationNesterovGreedyResearch(self,steps,ranking,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHeavyBallGreedyResearch(steps,ranking,gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationPolyakGreedyResearch(self,steps,gradTol,nMaxIter,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiPolyakResearch(steps,gradTol,nMaxIter,parallelism,nbThreads)
        
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime

    def perfomOptimisationNewResearch(self,gradTol,nMaxIter,w1,w2,wMaxIter,memSize=15,leven=1e-5,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiBFGSUnleashDaBeastResearch(gradTol,nMaxIter,w1,w2,wMaxIter,memSize,leven,parallelism,nbThreads)
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime
    def perfomOptimisationLM(self,gradTol,nMaxIter,w1,w2,wMaxIter,memSize=15,leven=1e-5,parallelism=False,nbThreads=-1):
        self.__Tri.perfomOptimPsiHessLM(gradTol,nMaxIter,w1,w2,wMaxIter,memSize,leven,parallelism,nbThreads)
        costf = np.array(self.__Tri.getCostFunc())
        gNorm = np.array(self.__Tri.getNorms())
        steps = np.array(self.__Tri.getSteps())
        times = np.array(self.__Tri.getTimes())
        hidd  = np.array(self.__Tri.getHidden())

        nbItt = self.__Tri.getnbItt()
        wallTime = self.__Tri.getWallTime()
        
        return costf,gNorm,steps,times,hidd,nbItt,wallTime