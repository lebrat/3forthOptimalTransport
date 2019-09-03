import sys
sys.path.insert(0, "../codeSemiDiscrete")
sys.path.insert(0, "../IntersectionGrille")
import numpy as np


def oracle(LagInstance,psi,w):
	"""
	@brief      Compute the cost function (scalar), the mass of each laguerre cell
				the gradient and its L2 norm and the number of hidden Laguerre cell
	
	@param      LagInstance  instance of LaguerreTri, the coordinates should be
							 set and the polyLine too.
	@param      psi          (np.array) The dual variable which drives the size of laguerre cell
							 at psi* the g(coord,psi*)=W_2^2(polyline,dirac(coord))
	@param      w            (np.array) the mass for each Dirac
	
	@return     cf,mass,Bar,grad,\|grad\|_2,nbHidden
	"""

	mass,Bar,cost = LagInstance.compute(psi)
	Cost = np.sum(cost) + np.dot(psi,w-mass)
	grad = w-mass
	numberHidden = np.sum(np.abs(mass)<1e-14)
	return Cost, mass, Bar, grad, numberHidden

def regularizedOracle(polyLine,x,psi,w,sigma=10.,offset=1e-9):	
	from scipy.ndimage.filters import gaussian_filter
	from OTInterface import OT, PostProcess,calcul_masse
	from optim import optimalTransportHessLM
	import interGridUtils
	box = [0,1,0,1]
	n = [2048,2048]	
	grid = interGridUtilsprint.InterGrid(box,n)
	density = grid.computeAllIntersections(polyLine)
	density = density.reshape(n[0],n[1],order='C')
	density = np.rot90(density,3)
	density = np.fliplr(density)
	img = gaussian_filter(density,sigma)+offset
	img/=calcul_masse(img,True)
	ot = OT(img,box,True)
	M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian=None)
	numberHidden = np.sum(np.abs(M)<1e-14)
	return Cost,M,Bar,g,numberHidden


class BFGS() :
	def __init__(self,nb_stock_max=8) :
		self.nb_stock_max=nb_stock_max
		self.Xfamily=[]
		self.gradfamily=[]
	def add_direction(self,s,d,gradfNew,gradf) :
		if len(self.Xfamily) == self.nb_stock_max :
			self.Xfamily.pop(0)
			self.gradfamily.pop(0)
		self.Xfamily.append(np.copy(s*d))
		self.gradfamily.append(np.copy(gradf-gradfNew))
	def find_direction(self, gradf) :
		if len(self.Xfamily)==0:
			return np.copy(gradf)
		else:
			r=np.copy(-gradf)
			stock=len(self.Xfamily)
			I=range(stock)
			rho=np.array([1./np.dot(self.Xfamily[i],self.gradfamily[i]) for i in I])
			alpha=np.zeros(stock)
			if stock > 0 :
				gamma=np.dot(self.Xfamily[-1],self.gradfamily[-1])/np.dot(self.gradfamily[-1],self.gradfamily[-1])
			else :
				gamma=1./np.linalg.norm(d)
			for i in reversed(I):
				alpha[i]=rho[i]*np.dot(self.Xfamily[i],r)
				r=r-alpha[i]*self.gradfamily[i]
			r=gamma*r
			for j in I:
				beta=rho[j]*np.dot(self.gradfamily[j],r)
				r+=(alpha[j]-beta)*self.Xfamily[j]
			return -r        

if __name__ == '__main__':
	def Rosenbrock_oracle(X,a=1.,b=100.):
		f=(a-X[0])**2+b*(X[1]-X[0]**2)**2
		gradf=np.array([ -2*a+(2-4*b*X[1])*X[0]+4*b*X[0]**3, 2*b*(X[1]-X[0]**2)])
		Hessf=np.array([[2-4*b*X[1]+12*b*X[0]**2,-4*b*X[0]]
						,[-4*b*X[0],2*b]])
		return -f,-gradf, -Hessf
	LBFGS=BFGS(10)
	x0=np.array([0.,0.])
	s=1.
	f,gradf,Hessf=Rosenbrock_oracle(x0)
	i=0
	while np.linalg.norm(gradf)>1.e-5 and i<2000:
		i+=1
		d=LBFGS.find_direction(gradf)
		f,gradfNew,Hessf=Rosenbrock_oracle(x0+s*d)
		LBFGS.add_direction(s,d,gradfNew,gradf)
		x0=x0+s*d
		print('x=(%1.3e,%1.3e) -- i:%3i -- Cost:%-02.6e -- ||g||:%1.6e -- alpha:%1.2e-- montee ?:%1.2e' %(x0[0],x0[1],i,f,np.linalg.norm(gradfNew),s,np.dot(d,gradf)))
		gradf=gradfNew
