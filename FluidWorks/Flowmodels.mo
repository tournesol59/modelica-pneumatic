within FluidWorks;
package Flowmodels
//replaceable package SimpleMedium=FluidWorks.SimpleAir.SimpleAir;

  connector FluidPort //"Generic fluid port"
     replaceable package Medium=SimpleMedium;
     flow Real m_flow "Flow into connector";
     Real p "Thermodynamic pressure at the connector";
     stream Real h_outflow "Specific enthalpy of outgoing fluid";
  end FluidPort;

  model TwoPorts   
      replaceable package Medium=SimpleMedium;
      FluidPort port_a(redeclare package Medium=Medium);
      FluidPort port_b(redeclare package Medium=Medium);
  end TwoPorts;

  model ControlVolume
     extends TwoPorts;
     parameter Real Volume=0.005;
     Real p, h;
     equation
        Volume/287/288*der(p) = port_a.m_flow + port_b.m_flow;
        Volume*der(p)+der(h)=port_a.m_flow*actualStream(port_a.h_outflow) + port_b.m_flow*actualStream(port_b.h_outflow); //"volume*der(p)+der(h)"
        port_a.p=p;
        port_b.p=p;
        port_a.h_outflow=h;
        port_b.h_outflow=h;
  end ControlVolume;

  model Flowmodel
    extends TwoPorts;
    parameter Real Area=0.001;
    Real rho; // ha, hb;
    equation
       port_a.m_flow=Area*sqrt(2*rho*abs(port_a.p-port_b.p));
       rho=(port_a.p+port_b.p)/2/287/288;
       //ha=inStream(port_a.h_outflow);
       //hb=inStream(port_b.h_outflow);
       port_a.m_flow+port_b.m_flow=0;
       port_b.h_outflow=inStream(port_a.h_outflow);
       port_a.h_outflow=inStream(port_b.h_outflow);
  end Flowmodel;
end Flowmodels;
