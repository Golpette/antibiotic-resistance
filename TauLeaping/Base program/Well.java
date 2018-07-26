import java.util.HashMap;

public class Well{




    // Ciprofloxacin concentration
    public double cipro;
    // Population number of each strain
    public HashMap<String, Integer> abundances ;

    // To hold changes in strain numbers due to (1) birth-death-mutation 
    // and (2) migration
    public HashMap<String, Integer> reproduction ;   // MESSY: this holds absolute number of events
    public HashMap<String, Double> migration ;       //        this holds expected events
  




    public Well(double c){  
        cipro = c;
        abundances = new HashMap<String, Integer>();

        reproduction = new HashMap<String, Integer>();
        migration = new HashMap<String, Double>();
    }

    




}
