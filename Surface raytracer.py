c=-4.1
def surface(x,y,z):
	#points on this surface yield 0 for the following expression
	return z+1/(((x*x+y*y)**2+c)**.5)#1/(x*x+y*y)+z*z-1

def norm(x,y,z):
	#gradient of the surface expression is normal to the surface
	Dx=-2*x*(x*x+y*y)/((x*x+y*y)**2+c)**1.5#-2*x/(x*x+y*y)**2
	Dy=-2*y*(x*x+y*y)/((x*x+y*y)**2+c)**1.5#-2*y/(x*x+y*y)**2
	Dz=1
	return Dx,Dy,Dz

def pattern(x,y,z):
	#color every point on the surface like so
	#return ((x>0),(y>0),(z>0))
	if (x%2-1)*(y%2-1)>0: #checkerboard pattern
		return (1,0,0) if z>=0 else (1,1,0)
	else:
		return (0,1,1) if z>=0 else (0,0,1)

#-------- The below functions work on any surface --------

def cross(v1,v2):
	#cross product of two vectors
	a,b,c=v1
	d,e,f=v2
	A=b*f-c*e
	B=c*d-a*f
	C=a*e-b*d
	return A,B,C

def unit(v1):
	#unit vector
	x,y,z=v1
	r=(x*x+y*y+z*z)**.5
	return x/r,y/r,z/r

def step(point,route,length):
	#approximate Euler integration
	x,y,z=point
	grad=norm(x,y,z)
	newPoint=x+route[0]*length,y+route[1]*length,z+route[2]*length
	newGrad=norm(newPoint[0],newPoint[1],newPoint[2])
	newRoute=unit(cross(cross(grad,route),newGrad))
	return newPoint,newRoute,length

def snap(point,error,gamma):
	#bring point back to the surface
	x,y,z=point
	while abs(surface(x,y,z))>error:
		fixdir=unit(norm(x,y,z))
		fixlen=-gamma*surface(x,y,z)
		x+=fixdir[0]*fixlen
		y+=fixdir[1]*fixlen
		z+=fixdir[2]*fixlen
	return x,y,z

#-------- Let's actually render some stuff! --------
import canvas
from math import cos,pi

def gross(N,A,B,C,D):
	a,b,c=N
	sqrt=abs((-D*D+B*B+A*A)*c*c-2*(B*C*b+A*C*a)*c+(-D*D+C*C+A*A)*b*b-2*A*B*a*b+(-D*D+B*B+C*C)*a*a)**.5
	denom=(B*B+A*A)*c*c-2*(B*C*b+A*C*a)*c+(C*C+A*A)*b*b-2*A*B*a*b+(C*C+B*B)*a*a
	x=(b*(-C*sqrt-B*D*a)+c*(B*sqrt-C*D*a)+A*D*c*c+A*D*b*b)/denom
	y=(c*(-A*sqrt-C*D*b)+C*a*sqrt+B*D*c*c-A*D*a*b+B*D*a*a)/denom
	z=-(-A*b*sqrt+B*a*sqrt+(B*D*b+A*D*a)*c-C*D*b*b-C*D*a*a)/denom
	'''
	x=(b*(C*sqrt-B*D*a)+c*(B*sqrt-C*D*a)+A*D*c*c+A*D*b*b)/denom
	y=(c*(A*sqrt-C*D*b)-C*a*sqrt+B*D*c*c-A*D*a*b+B*D*a*a)/denom
	z=-(A*b*sqrt-B*a*sqrt+(B*D*b+A*D*a)*c-C*D*b*b-C*D*a*a)/denom
	'''
	return x,y,z

pos=-2,0,-1/(4+c)
N=norm(pos[0],pos[1],pos[2])
angle=0
A1,B1,C1=unit(cross(cross(N,(1,0,0)),N)) #3D x_hat projected onto local plane
#print('x_hat: ',A1,B1,C1)
D1=cos(angle)
A2,B2,C2=gross(N,A1,B1,C1,D1) #3D rotated x_hat
#print('x_hat_rot: ',A2,B2,C2)
num_steps=100

x_scale,y_scale=10,10
origin=3,5

size=240,240
canvas.set_size(size[0],size[1])

for x_px in range(size[0]):
	for y_px in range(size[1]):
		P=pos
		x_im=x_px/size[0]*x_scale-origin[0]
		y_im=y_px/size[1]*y_scale-origin[1]
		#print(x_im,y_im)
		if x_im!=0 or y_im!=0:
			ds=(x_im**2+y_im**2)**.5/num_steps
			D2=x_im/(x_im**2+y_im**2)**.5
			dir=gross(N,A2,B2,C2,D2)
			#print(dir)
			for _ in range(num_steps):
				#x0=(P[0]+origin[0])/x_scale*size[0]
				#y0=(P[1]+origin[1])/y_scale*size[1]
				#canvas.fill_pixel(x0,y0)
				P,dir,ds=step(P,dir,ds)
			#input()
			r,g,b=pattern(P[0],P[1],P[2])
			canvas.set_fill_color(r,g,b)
			canvas.fill_pixel(x_px,y_px)

canvas.save_png('GravityLens.png')
