within FluidWorks;
  package PassiveComponent extends GenericComponent;

//  import SimpleAir;
  import GenericComponent;
  replaceable package SimpleMedium=FluidWorks.SimpleAir.SimpleAir;

	redeclare model extends ComponentModel1Demo(Temperature_Ext=293) //"note on change la temperature exterieure juste pour l'example, il faudrait rajouter des lois dynamiques etc"
            	
        end ComponentModel1Demo;

        redeclare model extends ComponentModel2Demo(Temp_Ext=293) //"on change la temperature exterieure, voir plus haut ComponentModel1Demo"

        end ComponentModel2Demo;
  end PassiveComponent;

