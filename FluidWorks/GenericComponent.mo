within FluidWorks;
  package GenericComponent
//     import SimpleAir;
//   Modelica.Fluid.SimpleAir.SimpleAir SimpleMedium;  //the instance of the package in MSL dir
//   replaceable package SimpleMedium=FluidWorks.SimpleAir.SimpleAir;
  
   connector myPort
       replaceable package Medium=SimpleMedium;  
       Real p          "The unique variable (pressure, Pa) of this Fluid Port" ;
       Real m_flow     "Flow of matter; m_flow > 0 if flow into component"; 
       stream Real h_outflow;
   end myPort;

   partial model partialTwoPorts
       parameter Boolean flowreverse = true;
       replaceable package Medium=SimpleMedium;
       myPort port_a(redeclare package Medium = Medium);
       myPort port_b(redeclare package Medium = Medium);
   end partialTwoPorts;

   replaceable partial model ComponentModel1Demo
       extends partialTwoPorts;
       inner Real pressure_diff;
       inner Real m_flow;
       inner Real Temperature;
       inner Real density; 
       parameter Real Temperature_Ext=273;
       inner Real Energy_Transfer;
       Real Energy_mix_in;
       ComponentModel2Demo press_model(redeclare package Medium=Medium);
       equation

          port_a.m_flow+port_b.m_flow=0;
          m_flow=port_a.m_flow;
          pressure_diff=port_a.p-port_b.p;	
        // below equation should be corrected with the inStream operator
	//  port_a.m_flow*(Energy_mix_in)+port_b.m_flow*(port_b.h_outflow)+
	  port_a.m_flow*(if m_flow>0 then Energy_mix_in else port_a.h_outflow) + port_b.m_flow*(if m_flow<0 then Energy_mix_in else port_b.h_outflow) - press_model.Energy_Transfer=0;
 
	  inStream(port_a.h_outflow)=port_b.h_outflow;
	  inStream(port_b.h_outflow)=port_a.h_outflow;
         // port_a.h_outflow=inStream(port_b.h_outflow); // assume Energy_Transfr==0
	 // port_b.h_outflow=inStream(port_a.h_outflow);
          Temperature=FluidWorks.SimpleAir.SimpleAir.SimpleInterfaces.IdealGasMedium.T_ph(p=port_a.p, h=Energy_mix_in);
          density=FluidWorks.SimpleAir.SimpleAir.SimpleInterfaces.IdealGasMedium.rho_pT(p=port_a.p, T=Temperature);
         
   end ComponentModel1Demo;
  
    replaceable partial model ComponentModel2Demo
        replaceable package Medium=SimpleMedium;
	outer Real pressure_diff;
        outer Real m_flow;
        outer Real Temperature;
        outer Real density;
	outer Real Energy_Transfer;
        parameter Real area=0.001; 
        parameter Real Temp_Ext=273;
	parameter Real h_const=0.0;
	equation
	  m_flow=area*sqrt(2*density*pressure_diff);  
          Energy_Transfer=if (pressure_diff>0) then (Temp_Ext-Temperature)*h_const else 0.0;
    end ComponentModel2Demo;

  end GenericComponent;

