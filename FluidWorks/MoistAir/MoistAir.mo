within MoistAir;
partial package PartialMedium "NOT COMPLETE NOW!, Generic medium interface"
   constant Integer nX "number of substances";     

   replaceable partial model BaseProperties
       Real X[nX]; //"     ...   "
   end BaseProperties;    

   replaceable partial function dynamicViscosity
      input  Real p;
      output Real eta;
   end dynamicViscosity;

end PartialMedium;

package MoistAir "Special type of medium"
    extends PartialMedium(nX = 2);

     redeclare model extends BaseProperties (T(stateSelect=StateSelect.prefer))
       // replaces BaseProperties by a new implementation and
       // extends from Baseproperties with modification
       // note, nX = 2 (!)
    equation
        X = {0, 1};
    end BaseProperties;

     redeclare function extends dynamicViscosity
        // replaces dynamicViscosity by a new implementation and
        // extends from dynamicViscosity
    algorithm
       eta := 2*p;
     end dynamicViscosity;

end MoistAir; 

