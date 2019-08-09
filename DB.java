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
	public static double[][][] nthVaribleCorrelationMatrix(double[][] par1,double[][] par2){
		
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
	
	}
/**
 * Creates a single matrix for correlation between matrices
 * @param par1
 * @return
 */
public static double[][] nthCorrelationMatrix(double[][] par1){
		
	
		
		double[][]matrix = new double[4][4];
		
		double[][] order = new double[par1.length][par1[0].length];
		double[][] data = new double[par1.length][par1[0].length];
		for(int i = 0; i < par1.length; i++) {
			order[i] = order(par1[i]);
			data[i] = data(order[i]);
			
		}
		
		
			for(int i =0; i < 4; i++) {
				for(int j =0; j < 4; j++) {
					
					double[][] bounds = new double[par1.length][data[0].length]; 
					for(int x= 0; x<  par1.length; x++) {
						//for(int y= 0; y<  par1.length; y++) {
							bounds[x][i] = data[x][i];
							bounds[x][j] =data[x][i+1];
							
						//}
					}
					
					
					
					matrix[i][j] = IntersectCountMutliD(par1, bounds);
					}
			
		}
		
		
		
		return matrix;
}
	
	
	
	
	
	
	
	
	
	public static int IntersectCount(double[] par1, double[] par2, double[] boundsPar1, double[] boundsPar2) {
		HashMap<Integer, Double> par1Map = timesInBoundsWithIndex(par1, boundsPar1);
		HashMap<Integer, Double> par2Map = timesInBoundsWithIndex(par2, boundsPar2);
		int count  = 0;
		for(int i = 0; i < par1Map.size(); i++) {
			if(par2Map.containsKey(par1Map.keySet().toArray()[i])) {
				count++;
			}
			
		}
		return count;
	}
	
	/**
	 * For two martices
	 * @param par1
	 * @param par2
	 * @param boundsPar1
	 * @param boundsPar2
	 * @return
	 */
	public static int IntersectCountMutliD(double[][] par1, double[][] bounds) {
		
		int count  = 0;
		ArrayList<HashMap<Integer, Double>> par1Map = timesInBoundsWithIndexMutliD(par1, bounds);
		
		
		for(int i = 0; i < par1Map.size(); i++) {
			for(int j = i+1; j < par1Map.size() - 1; j++) {
				//System.out.println(i);
				for(int k = 0; k < par1Map.get(j).keySet().toArray().length; k++) {
					if(par1Map.get(i).containsKey(par1Map.get(j).keySet().toArray()[k])) {
						System.out.println(count);
						count++;
					}
				}
			}
			
		}
		
		return count;
	}
	/**
	 * For two martices
	 * @param par1
	 * @param bounds
	 * @return
	 */
	
	public static ArrayList<HashMap<Integer, Double>> timesInBoundsWithIndexMutliD(double[][] par1,double[][] bounds) {
		ArrayList<HashMap<Integer, Double>> out = new ArrayList<HashMap<Integer, Double>>();
		
		HashMap<Integer, Double> indexToVal = new HashMap<Integer, Double>();
		
		for(int i = 0; i < par1.length; i++) {
			for(int j = 0; j < par1[i].length; j++) {
				
				boolean p1 = par1[i][j] >= bounds[i][0];
				
				if(p1 && par1[i][j] <= bounds[i][1]) {
					
					indexToVal.put(i, par1[i][j]);
				}
			}
			out.add(indexToVal);// makes sure values are correctly copied
			indexToVal.clear();
		}
		
		return out;
	}	
	
	
	/**
	 * For two vectors
	 * @param par1
	 * @param bounds
	 * @return
	 */
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
	
	public static double[][] multiplicar(double[][] A, double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        double[][] C = new double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
	/**
	 * Does matrix subtaction of probilities in order to find true bayesian correlation
	 * @param A data set
	 * @param B control set
	 * @return A - B
	 */
	public static double[][] matrixSetSubtraction(double[][] A, double[][] B){
		double[][] C = new double[A.length][A[0].length];
		for(int i = 0; i < A.length; i++) {
			for(int j = 0; j < A[0].length; j++) {
				C[i][j] = A[i][j] - B[i][j];
				
			}
		}
		return C;
	}
	
	/**
	 * Gets the trace of a matrix
	 * @param par1
	 * @return
	 */
	public static double trace(double[][] par1) {
		double tr = 0;
		for(int i = 0; i < par1.length; i++) {
			tr += par1[i][i];
		}
		return tr;
		
	}
	/**
	 * gets correlation based on trace
	 * @param par1 matrix
	 * @return correlation
	 */
	public static double traceCorrelation(double[][] par1) {
		
		double sum =0.0;
		for(int i = 0; i < par1.length; i++) {
			for(int j = 0; j< par1[j].length; j++) {
				sum+= par1[i][j];
			}
		}
		
		return trace(par1)/sum;
		
	}
	
	/**
	 * get the bayesian correlation for each quartile
	 * @param par1 the matrix
	 * @return 
	 */
	public static double[] BayesianCorrelationQuartiles(double[][] par1) {
		double[] ans = new double[par1.length];
		
		for(int i =0; i < par1.length; i++) {
			double diag = par1[i][i];
			double pb = 0;
			for(int j = 0; j < i+1; j++) {
				pb += par1[j][i];
			}
			
			ans[i] = diag / pb;
		}
		
		
		return ans;
		
		
	}
	
	/**
	 *  gets the bayesian correlation for the line
	 * @param par1 the matrix
	 * @return
	 */
	
	public static double BayesianCorrelationLine(double[][] par1) {
		double bottom = 0;
		
		for(int i =0; i < par1.length; i++) {
			
			double pb = 0;
			for(int j = 0; j < i+1; j++) {
				bottom += par1[j][i];
			}
			
			
		}
		
		
		return trace(par1) / bottom;
		
		
	}
	
}
