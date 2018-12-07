import java.util.ArrayList;
import java.util.HashMap;

public class DB {

	/**
	 * creates a matrix based on quartiles
	 * @param par1 A or Y, unordered/functional input
	 * @param par2 B or X, unordered/functional input
	 * @return the matrix
	 */
	public static double[][] correlationMatrix(double[] par1, double par2[]){
		
		
		double[] orderX = order(par2);
		double[] orderY = order(par1);
		
		double [] dataX = data(orderX);
		double[] dataY = data(orderY);
		
		
		double[][]matrix = new double[4][4];
		
		
		
			for(int i =0; i < 4; i++) {
				for(int j =0; j < 4; j++) {
					
					matrix[i][j] = IntersectCount(par1, par2, new double[]{dataY[i], dataY[i+1]}, new double[]{dataX[j], dataX[j+1]});
					}
			}
			
		
		
		return matrix;
	}
	
	/**
	 * creates a correlation matrix based on matrices created by using word2vec
	 * @param par1 the
	 * 
	 * @return the matrix new matrix
	 */
	public static double[][][] nthCorrelationMatrix(double[][] par1,double[][] par2){
		
		double[][][]matrix = new double[par1.length][4][4];
		
		for(int e =0; e < par1.length; e++) {
			double[] row = new double[par1[e].length];
			double[] row2 = new double[par2[e].length];
			for(int x =0; x < par1[e].length; x++) {
				row[x] = par1[x][e];
				row2[x] = par2[x][e];
			}
			matrix[e] = correlationMatrix(row, row2);
		}
		
		return matrix;
		
		/**

		
		double[][]intersect = new double[4][4];
		for(int k = 0; k < par1.size(); k++) {
			
			for(int i =0; i < 4; i++) {
				for(int j =0; j < 4; j++) {
					
					double[] orderX = order(par2.get(k));
					double[] orderY = order(par1.get(k));
					
					double [] dataX = data(orderX);
					double[] dataY = data(orderY);
					
					intersect[i][j] = IntersectCount(par1.get(k), par2.get(k), new double[]{dataY[i], dataY[i+1]}, new double[]{dataX[j], dataX[j+1]});
					
					
					}
			}
			matrix.add(intersect);
		}
		
		
		return matrix;
		**/
	}
	
	
	
	
	
	
	
	
	
	
	
	public static int IntersectCount(double[] par1, double[] par2, double[] boundsPar1, double[] boundsPar2) {
		HashMap<Integer, Double> par1Map = timesInBoundsWithIndex(par1, boundsPar1);
		HashMap<Integer, Double> par2Map = timesInBoundsWithIndex(par2, boundsPar2);
		int count  =0 ;
		for(int i = 0; i < par1Map.size(); i++) {
			if(par2Map.containsKey(par1Map.keySet().toArray()[i])) {
				count++;
			}
			
		}
		return count;
	}
	
	
	public static HashMap<Integer, Double> timesInBoundsWithIndex(double[] par1,double[] bounds) {
		
		HashMap<Integer, Double> indexToVal = new HashMap<Integer, Double>();
		
		for(int i = 0; i < par1.length; i++) {
			if(par1[i] >= bounds[0] &&  par1[i] <= bounds[1]) {
				indexToVal.put(i, par1[i]);
			}
		}
		return indexToVal;
	}
	
	
	/**
	 * 
	 * @param par1 order par1
	 * @return 0 = min value, 1 = lower quartile, 2 = median, 3= upper quartile, 4 = max value
	 */
	public static double[] data(double[] par1) {
	double[] out = {lowestVal(par1), lowerQuartile(par1), median(par1),upperQuartile(par1), highestVal(par1)};
		
		
		return out;
	}
	
	
	
	/**
	 * order the input first
	 * @param par1 the input
	 * @return the Upper Quartile
	 */
	public static double upperQuartile(double[] par1) {
		ArrayList<Double> par2 = new ArrayList<Double> ();
		for(double i : par1) {
			par2.add(i);
		}
		ArrayList<Double> par3 = new ArrayList<Double> ();
		
		
		if(par1.length % 2 == 0) {
			for(int i = (int) (par1.length*.5); i  < par1.length; i++) {
				par3.add(par1[i]);
			}
		}
		else {
			for(int i = (int) ((1 + par1.length) * .5 -1); i  < par1.length; i++) {
				par3.add(par1[i]);
			}
			 
		}
		
		
		double[] quart = new double[par3.size()];
		for(int i = 0; i < par3.size(); i++) {
		quart[i] = par3.get(i);
		
		}
		
		
		return median(quart);
	}
	
	/**
	 * order the input first
	 * @param par1 the input
	 * @return the Lower Quartile
	 */
	public static double lowerQuartile(double[] par1) {
		ArrayList<Double> par2 = new ArrayList<Double> ();
		for(double i : par1) {
			par2.add(i);
		}
		ArrayList<Double> par3 = new ArrayList<Double> ();
		
		
		if(par1.length % 2 == 0) {
			for(int i = 0; i  < (int) (par1.length*.5); i++) {
				par3.add(par1[i]);
				//System.out.println(par1[i]);
			}
		}
		else {
			for(int i = 0; i  < (int) ((1 + par1.length) * .5 ); i++) {
				par3.add(par1[i]);
			}
			 
		}
		
		
		double[] quart = new double[par3.size()];
		for(int i = 0; i < par3.size(); i++) {
		quart[i] = par3.get(i);
		
		}
		
		
		return median(quart);
	}
	
	
	
	
	
	/**
	 * order the input first
	 * @param par1 the input
	 * @return the median
	 */
	public static double median(double[] par1) {
		
		if(par1.length % 2 == 0) {
			return (par1[par1.length/2] + par1[par1.length/2 - 1]) *.5;
		}
		else {
			return par1[(int) ((1 + par1.length) * .5) - 1];
		}
		
	}
	
	
	
	
	/**
	 * orders the set from smallest to largest
	 * @param par1
	 * @return
	 */
	public static double[] order(double[] par1) {
		double[] clone = par1.clone();
		double[] ans = new double[par1.length];
		
		for(int i = 0; i < ans.length; i++) {
			ans[i] = replacelowestVal(clone);
		}
		
		return ans;
		
		
		
	}
	
	public static double lowestVal(double[] par1) {
		double low = highestVal(par1);
		double[] clone = par1.clone();
		for(double i : clone) {
			
			if(i<low)
				low=i;
		}
		return low;
	}
	
	public static double[] lowestValArrayList(ArrayList<Double>par1) {
		double low = Integer.MAX_VALUE;
		int index = 0;
		for(int i =0; i < par1.size(); i++) {
			
			if(par1.get(i)<low) {
				low=par1.get(i);
				index = i;
			}
		}
		double[] re ={low, index};
		return re;
	}
	
	
	
	/**
	 * 
	 * @param par1 the array
	 * @return
	 */
	public static double replacelowestVal(double[] par1) {
		
		ArrayList<Double> par2 = new ArrayList<Double> ();
		for(double i : par1) {
			par2.add(i);
		}
		
		double low = highestVal(par1);
		
		
		//for(int i =0; i < par1.length; i++) {
			
			
		low=lowestValArrayList(par2)[0];
		
		
		int index = (int)(lowestValArrayList(par2)[1]);
		par2.set(index, highestVal(par1));
		for(int i =0; i < par2.size(); i++) {
			par1[i] = par2.get(i);
		}
			
		//}
		return low;
	}
	
	
	
	public static double highestVal(double[] par1) {
		double high = 0;
		double[] clone = par1.clone();
		for(int i =0; i < par1.length; i++) {
			
			if(par1[i] > high)
				high= clone[i];
		}
		return high;
	}
	
	
}
