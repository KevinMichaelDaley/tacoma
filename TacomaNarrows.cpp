#include <cmath>
#include <array>

double sgn(double x){
    if(x==0) return 0;
    if(x>0) return 1;
    return -1;
}
struct TacomaParameters{

    double bridge_damping;
    double support_damping;
    double support_intrinsic_frequency;
    double bridge_intrinsic_frequency;
    double bridge_forcing_frequency;
    double bridge_forcing_amplitude;
    double vortex_shedding_strength;
    double support_to_bridge_coupling;
};
double characteristic(TacomaParameters params,double x){//characteristic equation to find the normal modes of the linear coupled system
    double lambda=params.support_damping;
    double w=params.support_intrinsic_frequency;
    double mu=params.support_to_bridge_coupling;
    double A=params.bridge_forcing_amplitude;
    double h=params.bridge_damping;
    double Omega=params.bridge_intrinsic_frequency;
    
    return (x*x+lambda*x+w*w)*(x*x+h*x+Omega*Omega)-2*mu*mu*x*x*x*x;
}
double d_characteristic(TacomaParameters params,double x){//analytic derivative of the above for use with newton's method
    double lambda=params.support_damping;
    double w=params.support_intrinsic_frequency;
    double mu=params.support_to_bridge_coupling;
    double A=params.bridge_forcing_amplitude;
    double h=params.bridge_damping;
    double Omega=params.bridge_intrinsic_frequency;
    
    double derivative_product=   (2*x+lambda) * (x*x+h*x+Omega*Omega)//x'y
                                +(x*x+lambda*x+w*w) * (2*x+h);//xy'
    return derivative_product-8*mu*mu   *x*x*x;//power rule
}
double find_zero(TacomaParameters params){
    double b=params.support_intrinsic_frequency;
    double a=params.bridge_intrinsic_frequency;
    if(b>a){
	    std::swap(b,a);
    }
    a+=0.25;
    b-=0.25;
    if(b<0){
	    b=0;
    }
    //bisection method as a first pass
    for(int n=0; n<20; ++n){
        if(std::abs(b-a)<0.01){
            break;
        }
        double fa=characteristic(params,a);
        double fb=characteristic(params,b);
        double fc=characteristic(params, (b+a)/2.0);
        if(std::abs(fc)<0.0001){
            break;
        }
        if(sgn(fa)==sgn(fc)){
            a=(a+b)/2.0;
        } else{
            b=(b+a)/2.0;
        }
        
    }
    double x=(b+a)/2.0;
    double threshold=0.01;
    //newton's method with early exit on success
    for(int i=0; i<100; ++i){
        double f2=characteristic(params,x);
        if(std::abs(f2)<threshold){
            return x;
        }
        double fp2=d_characteristic(params,x);
        x-=f2/fp2;
    }
    return (b+a)/2.0;
}
        
        
TacomaParameters make_default_params(double bridge_intrinsic_frequency, 
                               double forcing_frequency,
                                           double bridge_forcing_amplitude){
                    TacomaParameters p= {0.05,
                            0.1,
                            0.95,
                            bridge_intrinsic_frequency,
                            forcing_frequency,
                            bridge_forcing_amplitude,
                            0.05,
                            0.0
                    };
                    return p;
}

void support_system_derivatives(TacomaParameters params, 
                                std::array<double,2> support_vel,
                                std::array<double,2> support_pos,
                                double bridge_accel,
				double t,
                                std::array<double,2>* support_out_accel,
                                bool linear,
                                double mismatch,
                                double eig
                               )
{
    double lambda=params.support_damping;
    double w= params.support_intrinsic_frequency;
    double A=params.bridge_forcing_amplitude;
    double mu=params.support_to_bridge_coupling;
    double gamma=(!linear) ? params.vortex_shedding_strength : 0.0;
    for(int i=0; i<2; ++i){
        double dx=support_vel[i];
        double x=support_pos[i];
        double ddv=bridge_accel;
        w+=mismatch*(1.0-2.0*i)/2.0;
        ((*(support_out_accel)))[i]=-lambda*dx-w*w*x-mu*ddv+(gamma+A*sin(t*params.bridge_forcing_frequency))*(!linear?sgn(dx):1.0);
    }
}

double bridge_accel(TacomaParameters params,
                    std::array<double,2> support_vel,
                    std::array<double,2> support_pos,
                    double bridge_vel,
                    double bridge_pos,
                    double t,
                    bool linear,
                    double mismatch,
                    double eig
                   )
{
    double v=bridge_pos;
    double dv=bridge_vel;
    std::array<double,2> ddx;
    support_system_derivatives(params,support_vel,support_pos,0,t,&ddx, linear,mismatch,eig);
    double sum=ddx[0]+ddx[1];
    double A=params.bridge_forcing_amplitude;
    double h=params.bridge_damping;
    double Omega=params.bridge_intrinsic_frequency;
    double gamma=(!linear) ? params.vortex_shedding_strength : 0;
    double mu=params.support_to_bridge_coupling;
    double r=1.0;
    return ((-2*h*dv-Omega*Omega*v-sum)+(gamma+A*sin(t*params.bridge_forcing_frequency))*(!linear?sgn(dv):1.0))/(1-mu*mu*2*r);
}


void derivatives(TacomaParameters params, std::array<double,6> state_in, std::array<double,6>* derivs, double time, bool linear, double mismatch, double eig){
  std::array<double,2> support_vel={state_in[0],state_in[2]};
  std::array<double,2> support_pos={state_in[1],state_in[3]};
  double bridge_vel=state_in[4];
  double bridge_pos=state_in[5];
  std::array<double,2> support_new_accel={0,0};
	  double bridge_new_accel=bridge_accel(params,support_vel,support_pos,bridge_vel,bridge_pos,time,linear, mismatch,eig);
  support_system_derivatives(params, support_vel, support_pos, bridge_new_accel,time, &support_new_accel,linear, mismatch,eig);
  derivs[0]={support_new_accel[0],
          support_vel[0],
          support_new_accel[1],
          support_vel[1],
          bridge_new_accel,
          bridge_vel};
}
void integrate_once(TacomaParameters params, std::array<double,6>* state, double time1, double time2, bool linear=false, double mismatch=0, double eig=0){
    double h=time2-time1;
    std::array<double,6> derivs;
    std::array<double,6> s1, s2, s3;
    std::array<double,6> k1,k2,k3,k4;
    derivatives(params, *state, &derivs, time1,linear,mismatch,eig);
    for(int i=0; i<6; ++i){
        k1[i]=derivs[i]*h;
        s1[i]=(*(state))[i]+k1[i]/2;
    }
    derivatives(params, s1, &derivs, time1+h/2,linear,mismatch,eig);
    for(int i=0; i<6; ++i){
        k2[i]=derivs[i]*h;
        s2[i]=(*(state))[i]+k2[i]/2;
    }
    derivatives(params, s2, &derivs, time1+h/2,linear,mismatch,eig);
    for(int i=0; i<6; ++i){
        k3[i]=derivs[i]*h;
        s3[i]=(*(state))[i]+k3[i];
    }
    derivatives(params, s3, &derivs, time1+h,linear,mismatch,eig);
    for(int i=0; i<6; ++i){
        k4[i]=derivs[i]*h;
    }
    for(int i=0; i<6; ++i){
        (*(state))[i]+=(k1[i]+2*(k2[i]+k3[i])+k4[i])/6.0;
    }
}

#include <sstream>
#include <iostream>
#include <fstream>
int main(int argc, char** argv){
    int T=10000   ;
    double h=0.01;
    int N=1024;
    
    double max_amplitude=-10000000;
   
#ifndef TRACE_ONLY
    TacomaParameters peak_params;
        /*{
        std::ofstream ofs("eig.txt");
        
		   printf("linear sweep\n");
		  fflush(stdout);
        double delta=1;
        for(int i=0; i<N; ++i){
                    
                    for(int k=0; k<N; ++k){
                    int zero_crossing1=0;
                    int phase_diff=0;
                    int period=0;
                    int period_avg=0;
                    double maximum=-10000, minimum=100000;
                    double Omega=0.98;
                    double Beta=k/double(N)-0.5+0.9;
                    double Ay=12*i/double(N);
			    TacomaParameters params=make_default_params(Omega, Beta, Ay);
			    std::stringstream s;
			    s<<Omega<<"_"<<Ay<<"_"<<delta<<"_"<<0.0001;
			    
			    double eig=std::abs(find_zero(params));
			    std::array<double,6> state={delta,0.0,delta+0.25,0,0.0001,0};
			    std::array<double,6> state0={delta,0.0,delta+0.25,0,0.0001,0};
			    for(int t=0; t<T; ++t){
				double time0=t*h;
				double time1=(t+1)*h;
				
				integrate_once(params,&state,time0,time1,true,0,eig);
				
				if (t>T/2){
				    maximum=std::max(maximum,state[5]);
				    minimum=std::min(minimum,state[5]);
				    if(state0[1]<0 && state[1]>=0){
					if(zero_crossing1>0){
					    period=(t-zero_crossing1);
					}
					zero_crossing1=t;
				    }
				    if(state0[3]<0 && state[3]>=0){
					if(zero_crossing1>0){    
					
						    phase_diff=t-zero_crossing1;
						    period_avg=period;
				       }
					
				    }
				    for(int i=0; i<6; ++i){
					state0[i]=state[i];
				    }
				}
			    }
			    double amplitude=(maximum-minimum);
				    if(amplitude>max_amplitude){
					max_amplitude=amplitude;
					peak_params=params;
				    }
				    
				    
				    double tau=period_avg*h;
					    double pdiff=std::abs(fmod(double(phase_diff*h/tau),1.0));
					    
					    ofs<<Omega<<" "<<Ay<< " "<<" "<<Beta<<" "<<amplitude<<" "<<pdiff<<" "<<eig<<"\n";
					    ofs.flush();
					}
				    }
			    }*/
			   printf("nonlinear sweep\n");
			  fflush(stdout);
			    {
				std::string filename="sweep_deck.txt";
				std::ofstream ofs(filename);

				
				    double delta=0.001;
				    for(int i=0; i<N; ++i){
					for(int k=0; k<N; ++k){
					    int zero_crossing=0;
					    int Nor=0;
					    int phase_diff=0;
					    int period=0;
					    int period_avg=0;
					    double maximum=-10000, minimum=100000;
					    double Omega=1.0;
					    double Ay=12*i/double(N);
					    double Beta=k/double(N)-0.5+0.9;
					    TacomaParameters params=make_default_params(Omega, Beta, Ay);
				    std::stringstream s;
				    s<<Ay<<"_"<<delta<<"_"<<0.0001;
				    
				    std::array<double,6> state={-0.12,0.0,0.1,0,0.0,0};
				    std::array<double,6> state0={-0.12,0.0,0.1,0,0.0,0};
				    double orcos=0, orsin=0;
					double tmx=0, tmn=T;
					int crossings=0;
				    for(int t=0; t<T; ++t){
					double time0=t*h;
					double time1=(t+1)*h;
					
					integrate_once(params,&state,time0,time1);
					if(state[5]>0 && state0[5]<=0){
						tmn=std::min(tmn,t*h);
						tmx=std::max(tmx,t*h);
						crossings+=1;
					}
					if (t>T/2){
					    maximum=std::max(maximum,state[5]);
					    minimum=std::min(minimum,state[5]);
					    if(state0[1]<0 && state[1]>=0){
						if(zero_crossing>0){
							period=t-zero_crossing;
						}
						zero_crossing=t;
				    }
				    if(state0[3]<0 && state[3]>=0){
					if(period>0){
					    double diff=fabs(t-zero_crossing)*2*M_PI;
					    diff=fmod(diff,2*M_PI*period)/double(period);
					    phase_diff+=std::min(diff,2*M_PI-diff);
					    ++Nor;
					}

				    }
				}    
				
				    for(int i=0; i<6; ++i){
					state0[i]=state[i];
				    }
				
			    }
			    double amplitude=maximum-minimum;
			    double ord=phase_diff/(double)Nor;
			    double pdiff=ord;
			    ofs<<" "<<Ay<< " "<<" "<<Beta<<" "<<amplitude<<" "<<pdiff<<" "<<crossings/(tmx-tmn)<<"\n";
			    ofs.flush();
			    }
		       
		    }
		}
/*	    peak_params=make_default_params(0.99,1.0,9.2);
	    {
		   printf("basin sweep\n");
		  fflush(stdout);
		std::ofstream ofs_basin("basin.txt");
		for(int i=0; i<N; ++i){
		    for(int j=0; j<N; ++j){
			
			    double maximum=-100000,minimum=10000;
			    int zero_crossing1=0, period=0, period_avg=0;
			    int phase_diff=0;
			    double dx=i/double(N)-0.5;
			    double dy=j/double(N)-0.5;
			    double Omega=1.0;
			    double Ay=9.2;
			    double Beta=1.0;
			    std::array<double,6> state={dx,0.0,dy,0.0,0.0,0};
			    std::array<double,6> state0={dx,0.0,dy,0.0,0.0,0};
			    double orcos=0,orsin=0;
			    int Nor=0; 
			    int T2=T;
			    for(int t=0; t<T2; ++t){
				double time0=t*h;
				double time1=(t+1)*h;
				
				integrate_once(peak_params,&state,time0,time1);
				
				
				if (t>T2/2){
				    maximum=std::max(maximum,state[5]);
				    minimum=std::min(minimum,state[5]);
				    if(state0[1]<0 && state[1]>=0){
					    if(zero_crossing1>0){
						   
						    period=t-zero_crossing1;
					    }
					    
					zero_crossing1=t;
				    }
				    if(state0[3]<0 && state[3]>=0){
					if(period>0){
					    double diff=fabs(t-zero_crossing1)*2*M_PI;
					    diff=fmod(diff,2*M_PI*period)/double(period);
						    phase_diff+=std::min(diff,2*M_PI-diff);
						    ++Nor;
						}

					    }
					}    
					
					    for(int i=0; i<6; ++i){
						state0[i]=state[i];
					    }
					
				    }
				    double ord=phase_diff/(double)Nor;
				    double pdiff=ord;
				    
				    double amplitude=maximum-minimum;
				    ofs_basin<<dx<<" "<<dy<< " "<<" "<<Beta<<" "<<amplitude<<" "<<pdiff<<"\n";
				    ofs_basin.flush();
			    }
			}
		    }
		    {
		   printf("mismatch sweep\n");
		  fflush(stdout);
			   {
			    std::ofstream ofs("mismatch5.txt");
			
			    double max_amplitude=-10000000;
			    
			    double delta=0.1;
			    for(int i=0; i<N; ++i){
				    
				for(int k=0; k<N; ++k){
				    int zero_crossing1=0;
				    int phase_diff=0;
				    int period=0;
				    int period_avg=0;
				    double maximum=-10000, minimum=100000;  
				    double mismatch=i/double(N);
				    double Ay=12;
				    double Beta=k/double(N)*1-0.5+0.9;
				    double Omega=0.98;
				    TacomaParameters params=make_default_params(Omega, Beta, Ay);
			    
				    
				    
			    std::array<double,6> state={delta,0.0,delta+0.25,0,0.0001,0};
			    std::array<double,6> state0={delta,0.0,delta+0.25,0,0.0001,0};
				    double orcos=0,orsin=0;
				    int Nor=0; 
				    for(int t=0; t<T; ++t){
					double time0=t*h;
					double time1=(t+1)*h;
					
					integrate_once(params,&state,time0,time1,false,mismatch);
					
					
					if (t>T/2){
					    maximum=std::max(maximum,state[5]);
					    minimum=std::min(minimum,state[5]);
					    if(state0[1]<0 && state[1]>=0){
						    if(zero_crossing1>0){
							   
							    period=t-zero_crossing1;
						    }
						zero_crossing1=t;
					    }
					    if(state0[3]<0 && state[3]>=0){
						if(period>0){
						    double diff=fabs(t-zero_crossing1)*2*M_PI;
						   
					    		diff=fmod(diff,2*M_PI*period)/double(period);
						    phase_diff+=std::min(diff,2*M_PI-diff);
                                    ++Nor;
				}

                            }
                            
                        }
                            for(int i=0; i<6; ++i){
                                state0[i]=state[i];
                            }
                        
                    }
		    double ord=phase_diff/(double)Nor;
                    double amplitude=maximum-minimum;
                    
                    double pdiff=ord;
                    ofs<<" "<<mismatch<< " "<<" "<<Beta<<" "<<amplitude<<" "<<pdiff<<"\n";
                    ofs.flush();
                
            }
        }
    }
	}    */
   #else
		    { 
				    std::array<double,6> state={0.001*(rand()%20-10),0.0,0.001*(rand()%20-10),0,0.0,0};
				    double Ay=9.2;
				    double Beta=0.95;
				    double mismatch=0.0;
				    double Omega=0.98;
				    TacomaParameters params=make_default_params(Omega, Beta, Ay);
				    
				    fflush(stderr);
				    double orcos=0,orsin=0;
				     double At=0;
				     double h=0.01;
				     for(int t=0; t<100000; ++t){
					    	fprintf(stderr, "%lf %lf %lf %lf\n", t*h, state[1], state[3], state[5]);
					double time0=t*h;
					double time1=(t+1)*h;

					integrate_once(params,&state,time0,time1,false,mismatch);
					double theta=std::abs(state[3]-state[1]);
					At+=theta/10000;

				     }
				     printf("%f %f %f\n", Beta,mismatch, At);
				    fflush(stdout); 
		}	   

   #endif
}
            
                
                
            
            
        
