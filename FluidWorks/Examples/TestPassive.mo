within FluidWorks.Examples;
package TestPassive
  import FluidWorks.SimpleAir;
  import FluidWorks.GenericComponent;//"already imported by PassiveComponent"
  import FluidWorks.PassiveComponent;
 
    replaceable package mySimpleMedium=FluidWorks.SimpleAir.SimpleAir;  // the instance of the package in local dir

      model Generator
	 GenericComponent.myPort port_source(redeclare final package Medium=mySimpleMedium);
	 parameter Real Energy=500000;  //Joules
	 parameter Real p_upsm=1500000; //Pascal
         //Modelica.Blocks.Sources.Sine sinSource(freqHz=0.033);
	 equation
	     //port_source.m_flow; to be determined
	     port_source.p=p_upsm;
	     port_source.h_outflow=Energy;
      end Generator;


      model Ground
	  parameter Real p_amb=1000000; //Pascal
	  GenericComponent.myPort port_amb(redeclare final package Medium=mySimpleMedium);
      equation
	  port_amb.p=p_amb;
	  //port_amb.h_outflow=5000000; // NO!
	  
      end Ground;

      model useComponent
          Generator gen1;
	  PassiveComponent.ComponentModel1Demo tube(redeclare final package Medium=mySimpleMedium);
	  Ground gnd;
      equation
	  connect(gen1.port_source, tube.port_a);
	  connect(tube.port_b, gnd.port_amb); 
      end useComponent;

end TestPassive;
