within examples;
package TestMotor

  import Electric.Electricmotor;

  model useMotor
    Modelica.Blocks.Sources.Step step1(height=100);
//    Modelica.Blocks.Sources.Pulse pulse1(amplitude=100, width=50, period=0.1);
    Electric.Electricmotor.DCmotor motor1;
    Electric.Electricmotor.firstFilter filter1(moment_const=0.01);
    Real speed_rpm;
  equation
    connect(step1.y, motor1.command);
    connect(motor1.N_rpm, filter1.speed);
    connect(filter1.moment, motor1.resistive);
    speed_rpm=motor1.N_rpm;
  end useMotor;

end TestMotor;
