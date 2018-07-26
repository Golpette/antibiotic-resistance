import java.util.Arrays ;
import java.io.* ;
import java.util.HashMap;


// PROOF OF PRINCIPLE.
// MC optimization of diffusion parameters for each connecting channel. 
// Tau-leaping algorithm to simulate growth, mutation and migration  as in Bartek's experiments.

// NOTES:  (some things have changed from initial code)
//      - food is now explicitly present again. Without it, when diffusion params changed
//          per well we were increasing carrying capacity of system since any well below
//          carrying capacity was able to reproduce (i.e. infinite food). 
//          Now have 1-food in getGrowthRate().




public class TauLeap_MCopt_diff_FOOD_Channels{


    // ################## PARAMETERS ####################

    static final int num_optimisations = 200;

    static final int nwells = 24;  

    static final double timestep = 0.001 ;  // fix these 2 values for 4 days = 96 hours
    static final int numsteps = 40000 ;

    static final int K = 1000000 ;
    //static final int n_init = 10000;  //using K/2
    static int[] food = new int[nwells] ;

    static double[] D = new double[nwells-1] ;   //now each channel has specific rate (w: 1->2 = 2->1)

    static final double cipro_uni = 0;
    static double cipro_max = 0;
    //static final int plateauLength = 0;

    static final double bpMutation = 0;
    //static final double bpMutation = 0;
    static final double marRfact = 0;
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


    public static double getGrowthRate( double food_dens, double[] NK_array, double[] LB_rates ){
	// Method to estimate growth rate at specific N/K from Bartek's experimental data

        double NK = 1.0-food_dens;  // SUBSTITUTING pop density for food density

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




    //TODO : methods to calculate "distance" between simulation and experiment

    public static double getDist_timePoint(double[] simulation, Double[] experiment){
        //method to compare 2 sets of well densities at a given points in time  TODO:NOT IDEAL, requires time to match to 3 decimal places between data and output
        double dist = 0;
        if( experiment.length != nwells  || simulation.length != nwells ){ 
            System.out.println("Data being compared not consistent. Exiting");  System.exit(1);
        }

        for( int w=0; w<nwells; w++ ){
            dist += (experiment[w]-simulation[w])*(experiment[w]-simulation[w]);
        }

        //System.out.println( "dist = " + dist );
        return dist;
    }

















    public static void main(String args[])throws IOException{


        double BEST_DIST = 1E20;
        double[] BEST_PARAMS = new double[nwells];



        // Data output -------------
        PrintWriter out = new PrintWriter( new FileWriter("nwells.dat")  ) ;
        //PrintWriter out_antibiotic = new PrintWriter(new FileWriter("antibiotic.dat") );
        //PrintWriter out_finalWell = new PrintWriter( new FileWriter("finalWellPop.dat") );
        PrintWriter out22 = new PrintWriter( new FileWriter("clustering.dat")  ) ;


        // TODO Read in "experimental data" I am trying to reproduce; here I am using well densities vs time:
        //    Here: time series of well populations from
        HashMap<String, Double[]> EXP_DATA = Read.readExpData("TEST_DATA");




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



        //double reachedfinal_cutoff = 0.15 ; //i.e. if population is 15% of carrying capacity, say well is occuppied
        //int count_reachFinal = 0;


/**
        // TODO choose initial parameters
        double[] D_hold = new double[nwells-1];
        for( int aw=0; aw<nwells-1; aw++){
            //if( aw%2==0){  D_hold[aw] = 0.06 ; }
            //else{ D_hold[aw] = 0.02 ; }
            //D_hold[aw]=0.04;
            D_hold[aw]=0.08*Math.random();
            //D_hold[aw] = 0.08*Math.random() + 0.01; //(keep em same order of mag)
        }
**/

    
        ///// This set is "weird parameter set" with dist=2.66
        //double[] D_hold = {0.0254865423571532, 0.023816856517153942, 0.03846000255313947, 0.09002661081837596, 0.03417042838205398, 0.018145800879626434, 0.13667227085806666, 0.02355774574399946, 0.08003040782607855, 0.04656747311377386, 0.024461384150689554, 0.047267902339224184, 0.11427668671970717, 0.023452460665142477, 0.04312696434860049, 0.02137631075641685, 0.035600996628839865, 0.032448948898736396, 0.028652699635333737, 0.123971480291817, 0.10006632078104356, 0.014403649931679694, 0.01597771668669914, 0.0};


        ///// This is one with dist=147
        //double[] D_hold = {0.062239876971599745, 0.06965861604132467, 0.08761031823502588, 0.10705450949796064, 0.06691847054121121, 0.06981279012659128, 0.07074317144585392, 0.0383947360765151, 0.0792925209772357, 0.06956376610467258, 0.0342473983681339, 0.04812594179324564, 0.10648874473516258, 0.06351239321422394, 0.08402559533914045, 0.11725225688687804, 0.040602865496391576, 0.04466512947172466, 0.09735517480622802, 0.09324776639917733, 0.10020209981258865, 0.05373927962283667, 0.07287303122136551};



        double cum_dist_hold = 1E20;




	for( int exp=0; exp<num_optimisations; exp++ ){


        // TODO choose initial parameters
        double[] D_hold = new double[nwells-1];
        for( int aw=0; aw<nwells-1; aw++){
            //if( aw%2==0){  D_hold[aw] = 0.06 ; }
            //else{ D_hold[aw] = 0.02 ; }
            //D_hold[aw]=0.04;
            D_hold[aw]=0.2*Math.random();
            //D_hold[aw] = 0.08*Math.random() + 0.01; //(keep em same order of mag)
        }



	    // ----------------------- RESET EVERYTHING ---------------------
	    double t = 0;

	    Well[] wells = new Well[ nwells ];
	    for (int w=0; w<nwells; w++){
		    double cipro = cipro_uni + Math.pow( cipro_max, (1.0*w+1.0)/(1.0*nwells) );
		    wells[w] = new Well( cipro );
	    	//System.out.println( "well "+w+"  cipro="+cipro);
            //if( exp==0 ){  out_antibiotic.println(w + " " +cipro); }
            
            food[w] = K;  food[0]=K/2;
	    }
        //out_antibiotic.close();
	    // Initial inoculation of WT
	    wells[0].abundances.put("00000", K/2);     // NEED TO START THIS FROM SAME POINT AS EXPERIMENTAL DATA
    
        // cumulative distance between timepoints in since experiment
        double cum_dist = 0;

        // Set diffusion parameters from D_hold
        for( int www=0; www<nwells-1; www++){
            D[www] = D_hold[www];
        }
        
        // -------------------- RESET ------------------------------------



        // ------------------ TODO Modify parameter set -------------------
        for(int www=0; www<nwells-1; www++){
            //double random_shift = (D[www] * 0.5) * (Math.random() - 0.5) ;            // PROB: WE WILL ONLY SHIFT SMALL VALUES BY TINY AMOUNT WITH THIS
            //double random_shift = (Math.random() * 0.1 - 0.05) ;                    // Shift values by non-relative amount         ALLOWS NEGATIVE D VALUES   
            double random_shift = 0;
            D[www] += random_shift ;
        }
        // ----------------------------------------------------------------


        
        //System.out.println(Arrays.toString(food));
        


	    // EVOLVE SYSTEM for each experiment --------------------------------
	    for( int steps=0; steps<numsteps; steps++){


		// Get propensities for every genotype in every well at current time (i.e. before changing anything)
		for( int w=0; w<nwells; w++){

		    int well_pop = getWellPopulation( wells[w] );

//		    double well_dens = (double)well_pop/(double)K ;
// NOW USE FOOD DENSITY INSTEAD OF POPULATION DENSITY
            double food_dens = (double)food[w]/(double)K ;


		    // consider each strain separately
		    for( String key : wells[w].abundances.keySet() ){


			    int num_of_strain = wells[w].abundances.get( key ) ;

			    // BIRTH STUFF  (i.e. update wells[w]/ reproduction HashMap) 
			    //double g = getGrowthRate( well_dens, NKs, LB_growth_rates ) ;
                // *** USING FOOD DENSITY INSTEAD OF POPULATION DENSITY
			    double g = getGrowthRate( food_dens, NKs, LB_growth_rates ) ;

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
	    		wells[w].migration.put( key, num_of_strain * timestep ) ;    // HAVE REMOVED D[w] FROM HERE FOR CHANNEL METHOD *****
                           
             
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
				    to_left=0;  to_right = getPoissonRV( wells[w].migration.get( strain ) * D[w] ) ;   // HAVE PUT D[w] for move right, D[w-1] for left, HERE NOW *********
			    }                    
			    else if( w== nwells-1 ){
				    to_right=0;  to_left = getPoissonRV( wells[w].migration.get( strain ) * D[w-1] ) ;
			    }
			    else{
				    to_left = getPoissonRV( wells[w].migration.get( strain )  * D[w-1] ) ; 
				    to_right = getPoissonRV( wells[w].migration.get( strain ) * D[w]   ) ; 
			    }
			}
			else{
			    to_left = 0;   to_right=0;
			}

			// birth / death / mutation
			Integer currpop = wells[w].abundances.get(strain) ;
			if( currpop == null ){ currpop = 0; }       

			int add_fromBirth = wells[w].reproduction.get( strain ) ;

        	// Update food molecules
            food[w] -= add_fromBirth;

			int newpop = currpop + add_fromBirth - to_left - to_right ;                   
			if( newpop < 0 ){  newpop = 0 ; }   //don't bother getting new numbers
            
            // Update population
			wells[w].abundances.put( strain, newpop )  ; 




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




		// make nwells.dat file of well populations
		if(steps%100==0 && steps>0){

           //MAKE nwells.dat
		    out.print( String.format("%.3f", t)  );
		    for (int ww=0; ww<nwells; ww++){
			    //out.print( " "+getWellPopulation( wells[ww] )  );
			    out.print( " "+ (double)getWellPopulation( wells[ww] )/(double)K  );
		    }
		    out.println();


            double[] densities = new double[nwells];
            for(int g=0; g<nwells; g++){ 
                densities[g] = (double)getWellPopulation(wells[g])/(double)K ;  
            }


            // calculate distance between simulation & experiment here
            String ttt = String.format("%.3f", t) ;
            if( EXP_DATA.containsKey(ttt) ){
                double dist = getDist_timePoint( densities, EXP_DATA.get( ttt ) );
                cum_dist += dist;
            }


		}




	    }//end of numsteps



/**
        int finalWellPop = getWellPopulation( wells[nwells-1] );
        if( (double)finalWellPop/(double)K > reachedfinal_cutoff ){
            count_reachFinal += 1 ;
        }

	    for( int w=0; w<nwells; w++ ){
		    //System.out.println( "well "+w+" : "+ wells[w].abundances );
	    }     
**/


        // Print 3D clustering plot:
        double biggestD = 0;
        for( int azz=0; azz<D.length; azz++){ if(D[azz]>biggestD){ biggestD=D[azz]; } }
        double smallestD = 1E30;
        for( int azz=0; azz<D.length; azz++){ if(D[azz]<smallestD){ smallestD=D[azz]; }  }
        double sumSqrDiff = 0;
        for( int azz=0; azz<D.length; azz++){
            sumSqrDiff += (D[azz]-0.04)*(D[azz]-0.04);
        }
        
        out22.println( cum_dist + " " + (biggestD-smallestD) + " " + sumSqrDiff/D.length );
        
        




        // UPDATE BEST FOUND SOLUTION
        if( cum_dist < BEST_DIST )  { 
            BEST_DIST = cum_dist ; 
            for( int x=0; x<nwells-1; x++){
                BEST_PARAMS[x] = D[x];
            }            
        }



        //System.out.println("cum dist = " + cum_dist);
        System.out.println( cum_dist);


        // TODO: reject or accept parameters - if accepting, must update D_hold -------------------------
        if( cum_dist < cum_dist_hold ){ // THEN ACCEPT CHANGE

            //System.out.println("Accept");
            
            cum_dist_hold = cum_dist;
            for(int q=0; q<nwells-1; q++){
                D_hold[q] = D[q] ;
            }// otherwise leave D_hold unchanged and then we will use
             // these previous parameters as our next starting point        
        }
        //else{ System.out.println("Reject");  }


        //System.out.println( Arrays.toString(D_hold) );



	}//num optimisations




//    System.out.println("");
//    System.out.println("Best distance = "+ BEST_DIST );
//    System.out.println("Best params = "+ Arrays.toString(BEST_PARAMS) );

//    System.out.println("food = " + Arrays.toString(food) );

    //////System.out.println("Cipro max= "+cipro_max +" :  Num exps="+ num_experiments +", reached final well=" + count_reachFinal + "  P="+(double)count_reachFinal/(double)num_experiments);

    out.close();

    out22.close();


	System.exit(0);



    }//main




}//class






