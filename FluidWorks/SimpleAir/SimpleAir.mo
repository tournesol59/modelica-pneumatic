within SimpleAir;
   package SimpleAir

     package SimpleInterfaces
	record ThermodynamicState
		Real p;
		Real T;		
 	end ThermodynamicState;

	package IdealGasMedium


          constant Real HeatCapacity=1004.5;
          constant Real HeatCapacityIntern=717.5;

          function T_ph
             parameter Real T_init=293.16;
             input Real p;
             input Real h;
             Real h_init;
             output Real T;
             algorithm
                h_init:=HeatCapacity*T_init;
                T:=(h-h_init)/HeatCapacityIntern;
	        return;
           
          end T_ph;

          function rho_pT
             input Real p;
             input Real T;
             output Real rho;
             algorithm
                rho:=p/287/T;
	        return;
           
          end rho_pT;


          function h_pT
             parameter Real T_init=293.16;
             input Real p;
             input Real T;
             Real T_init;
             output Real h;
             algorithm
                h_init:=HeatCapacity*T_init;
                h:=h_init+(T-T_init)*HeatCapacityIntern;
	        return;
           
          end h_pT;

        end IdealGasMedium;

     end SimpleInterfaces;

   end SimpleAir;

