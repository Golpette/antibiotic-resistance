import java.util.HashMap;
import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Read{



    public static HashMap<String, Double> readSDMap( String filename )throws IOException{
        // Reading in a String-Double map of fitnesses/MICS for genotypes
        // FILE MUST HAVE 2 COLUMNS:  String genotype, Double fitness/MIC

        HashMap<String, Double> map = new HashMap<String, Double>();

        BufferedReader br = new BufferedReader(new FileReader(filename) );
        try{
            String line = br.readLine();

            while( line != null ){
            // Get components and add to HashMap
            Scanner scan = new Scanner(line);
            // File must have 2 columns
            String genotype = scan.next();
            double value = scan.nextDouble();

            map.put(genotype,value);

            line = br.readLine();

            }    
        } finally {
            br.close();
        }

        return map;
    }





    public static double[] readGrowthData( String filename, String demand )throws IOException{

        ArrayList<Double> ar = new ArrayList<Double>();

        BufferedReader br = new BufferedReader( new FileReader(filename) );
        try{
            // Read first column for "densities" and 2nd for "growth rates"
            String line = br.readLine();
            while( line!=null ){

                Scanner scan = new Scanner(line);
                if( demand.equals( "density" ) ){
                    ar.add( scan.nextDouble() );
                }
                else if( demand.equals( "growthrate" ) ){
                    //ignore 1st column
                    String ignore = scan.next();
                    ar.add( scan.nextDouble() );
                }
                else{
                    System.out.println( "WRONG DEMAND in readGrowthData()" );  System.exit(0);
                }
                
                line = br.readLine();
            }

        }finally {
            br.close();
        }

        double[] rt = new double[ ar.size() ];       
        int i=0;
        for( Double d : ar ){
            rt[i] = d;
            i += 1;
        }
        return rt ;
    }




    public static HashMap<String, Double[]> readExpData(String filename)throws IOException{
        // File must be "time" followed by well densities

        HashMap<String, Double[]> data = new HashMap<String, Double[]>();


        BufferedReader br = new BufferedReader(new FileReader(filename) );
        try{
            String line = br.readLine();
            while( line!=null ){

                ArrayList<Double> densities = new ArrayList<Double>();

                Scanner scan = new Scanner(line);
                String time = scan.next();
                

                while( scan.hasNext() ){
                    densities.add( scan.nextDouble() );
                    //System.out.println("xx");
                }


                //convert ArrayList to array
                Double[] rt = new Double[ densities.size() ];       
                int i=0;
                for( Double d : densities ){
                    rt[i] = d;
                    i += 1;
                }

                data.put( time, rt);
                
                line = br.readLine();
            }


        }finally{
            br.close();
        }

        return data ;

    }






}
