#include "generators.h"

EdwardsGenerator::EdwardsGenerator(Torsion t,unsigned int from,unsigned int b)
 : CurveGenerator(t,b)
 {
	switch (T)
	{
		case Z12:
			C = new EllipticCurve("0","0","0","-12");
			G = new RationalPoint("-2","-4");
			A = 1;
		break;
		case Z2xZ8:
			C = new EllipticCurve("0","0","0","-8");
			G = new RationalPoint("12","40");
			A = 1;
		break;
		case Z6:
			C = new EllipticCurve("0","-1/2304","0","-5/221184");
			G = new RationalPoint("1/192","1/4608");
			A = -1;
		break;
		case Z8:
			C = new EllipticCurve("0","0","0","48");
			G = new RationalPoint("4","-16");
			A = -1;
		break;
		case Z2xZ4:
			C = new EllipticCurve("0","0","0","-11664/25");
			G = new RationalPoint("36","-864/5");
			A = -1;
		break;
	};

	initialize(from);
 }
 
void EdwardsGenerator::generate_base_point(RationalPoint& P)
{	
	Qrac s = R.X,t = R.Y,x,y;
	
	if (T == Z12)
	{
		x  = 8*t*((s^2)+12);
		x /= (s-2)*(s+6)*((s^2)+12*s-12);
		
		y  = 4*s*((s^2)-12*s-12);
		y /= (s-2)*(s+6)*((s^2)-12);
		
		y  = -y;	
		
	}
	else if (T == Z2xZ8)
	{
		Qrac a = 1/((s-9)/(t+25)+1);		
		Qrac b = 2*a*(4*a+1)/(8*(a^2)-1);

		x  = (2*b-1)*(4*b-3)/(6*b-5);
		y  = (2*b-1)*((t^2)+50*t-2*(s^3)+27*(s^2)-104);
		y /= ((t+3*s-2)*(t+s+16));
		
	}
	else if (T == Z6)
	{
		Qrac sig = (1-96*s)/(96*s);
		
		Qrac r   = t/(s^2);
		
		Qrac al = (sig^2)-5;
		Qrac bt = 4*sig;

		x = 2*r*sig/((sig-1)*(sig+5)*(sig*sig+5));
		y = ((al^3)-(bt^3))/((al^3)+(bt^3));
	
	}
	else if (T == Z8)
	{
		Qrac u = 2*s/t;
		Qrac v = (2*(s^3)-(t^2))/(t^2);
		
		x = 2*(u^2);
		y = (4*(u^4)-1)/v;
		
	}
	else if (T == Z2xZ4)
	{
		x = "1/3";
		
		y = 625*(t^2);
		y = -y + 77760*t+1250*(s^3)+24300*(s^2)-3779136;
		y = ((30*s-25*t+1944)^2)/y;
	
	}
	else return;
			
	P.set(x,y);
}
