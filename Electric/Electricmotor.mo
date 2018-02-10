within Electric;

package Electricmotor
//approximation by a DC linear electrical motor
  connector RealInput = input Real;
  connector RealOutput = output Real;
  
   model DCmotor
     import SI=Modelica.SIunits;
//     SI.Voltage command=100;
     RealInput command;
//     SI.Torque resistive=1.0;
     RealInput resistive;
//     SI.Torque positive;
     Real positive;
//     SI.Angle angle_rad;
//     SI.AngularVelocity speed_rads; 
     Real angle_rad, speed_rads;
//     Output SI.AngularVelocity N_rpm;
     RealOutput N_rpm;
//     SI.Amp i_c; //current
     Real i_c;
     parameter SI.Inductance L_c = 0.001;
     parameter SI.Resistance R_c = 100;
     parameter SI.Inertia J_m = 10;
     parameter Real K_c=0.1;
     equation
         positive=K_c*i_c;
         J_m*der(speed_rads)=-0.0001*speed_rads + 
                         positive - resistive;
         
         command=R_c*i_c + L_c*der(i_c) + K_c*speed_rads; // electric equation
         der(angle_rad) = speed_rads;
         
         N_rpm=speed_rads/3.1415*30;
     initial equation
         speed_rads = 0.0;
         angle_rad = 0.0;
         i_c = 0.0;
   end DCmotor;

   model firstFilter
      import SI=Modelica.SIunits;
      RealInput speed;
      RealOutput moment;  // should be an output
//      parameter SI.Torque moment_end;
      parameter Real moment_const;
   equation
      der(moment)=-0.1*moment+10*moment_const*speed;
   initial equation
      moment=0.1;
   end firstFilter;

end Electricmotor;
