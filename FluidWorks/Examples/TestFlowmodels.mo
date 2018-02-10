within FluidWorks.Examples;
package TestFlowmodels 
  import FluidWorks.SimpleAir.SimpleAir;
  import FluidWorks.Flowmodels;
 
    replaceable package mySimpleMedium=FluidWorks.SimpleAir.SimpleAir;  //"the instance of the package in local dir" 

      model Generator
	 Flowmodels.FluidPort port_source(redeclare final package Medium=mySimpleMedium);
	 parameter Real Energy=500000;  //"Joules"
	 parameter Real p_upsm=1500000; //"Pascal"
         //Modelica.Blocks.Sources.Sine sinSource(freqHz=0.033); //"dynamic test"
	 equation
	     //port_source.m_flow=0.08; // "NO TRY WITH CVol Only"
	     port_source.p=p_upsm;  // "L UN OU L AUTRE"
	     port_source.h_outflow=Energy;
      end Generator;

      model Ground
	  parameter Real p_amb=1000000; //"Pascal"
	  Flowmodels.FluidPort port_amb(redeclare final package Medium=mySimpleMedium);
      equation
	  port_amb.p=p_amb;
	  //port_amb.h_outflow=5000000; //"NO!"
	  
      end Ground;

      model useComponent
          Generator gen1;
	  Flowmodels.ControlVolume tube(redeclare final package Medium=mySimpleMedium);
	  Flowmodels.Flowmodel valve(redeclare final package Medium=mySimpleMedium);
 	  Flowmodels.Flowmodel orifice(redeclare final package Medium=mySimpleMedium);
          Ground gnd;
      equation
	  connect(gen1.port_source, valve.port_a);
          connect(valve.port_b, tube.port_a);
          connect(tube.port_b, orifice.port_a);
	  connect(orifice.port_b, gnd.port_amb); 
      initial equation
          tube.p=101325;
          reservoir.p=101325;
          tube.h=1004*288;
          reservoir.h=1004*288;
      end useComponent;

end TestFlowmodels;
