import java.util.Arrays ;
import java.io.* ;
import java.util.HashMap;

// Tau-leaping algorithm to simulate growth, mutation and migration  as in Bartek's experiments.
//  


// Possible improvements:
//    (1) if number events expected is large, approximate Poisson distribution with normal distribution.
//    (2) 



public class TauLeapSystem{

    // ################## PARAMETERS ####################

    static final int num_experiments = 1;

    static final int nwells = 24;  

    static final double timestep = 0.001 ;  // fix these 2 values for 3 days = 72 hours
    static final int numsteps = 72000 ;

    static final int K = 10000000 ;
    static final int n_init = 10000;

    static final double D = 0.04 ;

    static final double cipro_uni = 0;
    static double cipro_max = 1000;
    //static final int plateauLength = 0;

    static final double bpMutation = 1E-8;
    //static final double bpMutation = 0;
    static final double marRfact = 10;
    // Mutation rate 'mu' depends on Cipro concentration, a, with form:  mu = bpMutation * [ 1 + ( X.a/MIC_WT )^Y  ]
    //    max_increase*bpMutation is upper limit on mutation probability
    //static double X = 0.0;   static final double Y = 1.0;  
    static double X = 15.0;   static final double Y = 1.5;  
    static final int max_increase = 200;                     

    // ##################################################




    public static int getPoissonRV( double lambda ){
	// Get number of events from Poisson distribution. Is very slow for large number of events.
        double l = Math.exp( -lambda ) ;
        double p = 1.0 ;  int k = 0 ;
        do{
            k++ ;
            p = p * Math.random() ;
        } while ( p > l ) ;
        return k-1 ;
    }


    public static int getBinomialRV(int num_trials, double prob) {
	// Get number of events from binomial distribution.
        int events = 0;
        for(int i = 0; i < num_trials; i++) {
            if(Math.random() < prob)
		events++;
        }
        return events;
    }


    public static double getGrowthRate( double NK, double[] NK_array, double[] LB_rates ){
	// Method to estimate growth rate at specific N/K from Bartek's experimental data
        double g = -9;
        int entry = 0;
        for( int i=0; i<NK_array.length; i++ ){
            if( NK_array[i] > NK ){
                entry = i;
                i = NK_array.length * 2 ;
            }
        }
        if( NK == 0 ){ g = 2; }
        else if( NK >= 1 ){ g = 0; }
        else{
            g = (LB_rates[entry] + LB_rates[entry-1]) / 2.0 ;
        }
        if( g < 0 ){ System.out.println("g<0: ERROR"); System.exit(0); }

        return g;

        // TURN OFF BIRTHS HERE
        //return 0;
    }



    
    public static int getWellPopulation( Well wl ){
	// Get total population of a well  -- THIS SHOULD BE IN WELL CLASS
        int totPop = 0;
        for( String key : wl.abundances.keySet() ){
            totPop += wl.abundances.get(key) ;
        }
        return totPop;
    }




    public static void implementMutations( String s, int num_reps, HashMap<String, Integer> repro, double MIC_WT, double antibiotic ){  //SHOULD IT BE WELL OBJECT HERE INSTEAD OF REPRO!?!?
	// modifies "reproduction" counter in Well object directly

        // Calculate mutation rate in this well
        double mutation_rate = bpMutation * ( 1 + Math.pow( X*antibiotic/MIC_WT, Y)  ) ;
        if( mutation_rate > bpMutation*max_increase ) { mutation_rate = bpMutation*max_increase ; }


        int counter = 0;
        for( int i=0; i < s.length(); i++){

            String strain = s ;
            char[] ca = strain.toCharArray();
            ca[i] = '1';
            String new_strain = String.valueOf( ca );


            if( new_strain.equals( strain ) ){
                //do nothing
            }
            else{
                // make mutants
                int num_muts;  // usually 0 or maybe 1
                if( i<3 ){
                    //num_muts = getPoissonRV( bpMutation * num_reps ) ;
                    num_muts = getPoissonRV( mutation_rate * num_reps ) ;
                }
                else{
                    // ##########   last 2, marR and acrR are more likely  ###########
                    //num_muts = getPoissonRV( marRfact * bpMutation * num_reps ) ;   
                    num_muts = getPoissonRV( marRfact * mutation_rate * num_reps ) ;                     
                }
                // count total number of mutants generated in this method
                counter = counter + num_muts;
                if( num_muts > 0 ){
                    repro.put( new_strain, num_muts );
                }
            }
        }

        // make sure we're not adding too many new bacteria
        int diff = num_reps - counter ;
        if( diff < 0 ){  diff = 0 ; }  //i.e. all bacteria birth events produced a mutant
        else{
            repro.put( s, diff ) ;
        }


    }









    public static void main(String args[])throws IOException{



        cipro_max = Double.parseDouble(args[0]);



        // Data output -------------
        PrintWriter out = new PrintWriter( new FileWriter("nwells.dat")  ) ;
        PrintWriter out_antibiotic = new PrintWriter(new FileWriter("antibiotic.dat") );
        PrintWriter out_finalWell = new PrintWriter( new FileWriter("finalWellPop.dat") );


        // Data input ------------------------
        // Read in growth data
        double[] NKs = Read.readGrowthData( "growthData_LB_nonMixed.dat", "density");
        double[] LB_growth_rates = Read.readGrowthData( "growthData_LB_nonMixed.dat", "growthrate");
        // Read in Map of fitnesses
        HashMap<String, Double> FITNESS = Read.readSDMap( "Fitnesses.dat" );
        // Read in Map of MICs 
        HashMap<String, Double> MICS = Read.readSDMap( "MICs.dat" );
        // Get wild type MIC
        double WT_MIC = MICS.get("00000");



        double reachedfinal_cutoff = 0.15 ; //i.e. if population is 15% of carrying capacity, say well is occuppied
        int count_reachFinal = 0;


	for( int exp=0; exp<num_experiments; exp++ ){



	    // RESET EVERYTHING ---------------------
	    double t = 0;

	    Well[] wells = new Well[ nwells ];
	    for (int w=0; w<nwells; w++){
		    double cipro = cipro_uni + Math.pow( cipro_max, (1.0*w+1.0)/(1.0*nwells) );
		    wells[w] = new Well( cipro );
	    	//System.out.println( "well "+w+"  cipro="+cipro);
            if( exp==0 ){  out_antibiotic.println(w + " " +cipro); }
	    }
        out_antibiotic.close();
	    // Initial inoculation of WT
	    wells[0].abundances.put("00000", n_init);
        
        




	    // EVOLVE SYSTEM for each experiment --------------------------------
	    for( int steps=0; steps<numsteps; steps++){


		// Get propensities for every genotype in every well at current time (i.e. before changing anything)
		for( int w=0; w<nwells; w++){

		    int well_pop = getWellPopulation( wells[w] );
		    double well_dens = (double)well_pop/(double)K ;



		    // consider each strain separately
		    for( String key : wells[w].abundances.keySet() ){


			    int num_of_strain = wells[w].abundances.get( key ) ;

			    // BIRTH STUFF  (i.e. update wells[w]/ reproduction HashMap) 
			    double g = getGrowthRate( well_dens, NKs, LB_growth_rates ) ;

			    //System.out.println( "g: "+g );

    			double fitness = FITNESS.get( key ) ;
	    		double birthPROP = num_of_strain * g * fitness ;
	    		double avg_birth_attempts = timestep * birthPROP ;                        


	    		// PROB OF DEATH upon_replication : (cipro/MIC)^2
	    		double probDeath = Math.pow( (wells[w].cipro/MICS.get(key)), 2) ;
	    		if( probDeath > 1 ){  probDeath = 1.0;  }


	    		//int rand_born = getPoissonRV( avg_birth_attempts * (1-probDeath) );
	    		int rand_born = getPoissonRV( avg_birth_attempts );                   
	    		///int rand_die  = getPoissonRV( avg_birth_attempts * probDeath     );        
	    		int rand_die = 2 * getBinomialRV( rand_born, probDeath );


	    		// Ensure populations don't become negative
	    		int net_births = rand_born - rand_die ;
	    		if( net_births + num_of_strain < 0 ){  net_births = -1 * num_of_strain ; }                

	    		// Store population increments
	    		if( net_births <= 0 ){
	    		    wells[w].reproduction.put( key, net_births ) ;
	    		}
	    		else{
	    		    //  MUTATIONS
	    		    implementMutations( key, net_births, wells[w].reproduction, WT_MIC, wells[w].cipro );
	    		}


	    		// MIGRATION STUFF  (i.e. update wells[w].migration HashMap)
	    		wells[w].migration.put( key, num_of_strain * D * timestep ) ;    
	    		// only putting avg here and will calculate 2 random numbers from it (LEFT AND RIGHT). This is stupid. //TODO
                           
             
		    }//end of strain loop

		}// end of wells loop








		//  Get migration propensities and update all strain abundances across wells ----------------------------------
		for( int w=0; w<nwells; w++ ){

		    for( String strain: wells[w].reproduction.keySet()  ){  // "reproduction" HashMap since may include new mutants!
			//   BUT some of these wont exist in "migration"


			// migration 
			int to_left;  int to_right;
			if( wells[w].migration.containsKey( strain ) ){

			    if( w==0 && nwells==1 ){
				to_left = 0;   to_right=0;
			    }
			    else if( w==0 ){
				to_left=0;  to_right = getPoissonRV( wells[w].migration.get( strain ) ) ;
			    }                    
			    else if( w== nwells-1 ){
				to_right=0;  to_left = getPoissonRV( wells[w].migration.get( strain ) ) ;
			    }
			    else{
				to_left = getPoissonRV( wells[w].migration.get( strain ) ) ; 
				to_right = getPoissonRV( wells[w].migration.get( strain ) ) ; 
			    }
			}
			else{
			    to_left = 0;   to_right=0;
			}


			// birth / death / mutation
			Integer currpop = wells[w].abundances.get(strain) ;
			if( currpop == null ){ currpop = 0; }       

			int add_fromBirth = wells[w].reproduction.get( strain ) ;

			int newpop = currpop + add_fromBirth - to_left - to_right ;                   
			if( newpop < 0 ){  newpop = 0 ; }   //don't bother getting new numbers
            
			wells[w].abundances.put( strain, newpop )  ; 
			//


			// Update neighoburing wells due to migration
			if( w > 0   &&   to_left > 0 ){
			    if( ! wells[w-1].abundances.containsKey( strain )) {
				wells[w-1].abundances.put(strain, to_left);
			    }
			    else {                            
				wells[w-1].abundances.put(strain, wells[w-1].abundances.get(strain)+to_left);
			    }
			}
			if( w < nwells-1   &&   to_right > 0){
			    if(! wells[w+1].abundances.containsKey( strain )) {
				wells[w+1].abundances.put(strain, to_right);
			    }
			    else {
				wells[w+1].abundances.put(strain, wells[w+1].abundances.get(strain)+to_right);
			    }
			}




		    }//end strains in reproduction HashMap loop 
		}//end nwells loop




		// clear all counters in each well
		for( int w=0; w<nwells; w++ ){
		    wells[w].reproduction.clear() ;
		    wells[w].migration.clear() ;
		}



		// increment time
		t += timestep ;


/**
		//long total_population = 0;
		// make nwells.dat file of well populations
		if(steps%10==0){
		    out.print( t );
		    for (int ww=0; ww<nwells; ww++){
			out.print( " "+getWellPopulation( wells[ww] )  );
			//total_population += getWellPopulation( wells[ww] ) ;
		    }
		    out.println();
		    //System.out.println(t + " " + total_population);
		}
**/


	    }//end of numsteps


        int finalWellPop = getWellPopulation( wells[nwells-1] );
        if( (double)finalWellPop/(double)K > reachedfinal_cutoff ){
            count_reachFinal += 1 ;
        }

	    for( int w=0; w<nwells; w++ ){
		    //System.out.println( "well "+w+" : "+ wells[w].abundances );
	    }     




	}//num experiments



    System.out.println("Cipro max= "+cipro_max +" :  Num exps="+ num_experiments +", reached final well=" + count_reachFinal + "  P="+(double)count_reachFinal/(double)num_experiments);

    out.close();

	System.exit(0);



    }//main




}//class






