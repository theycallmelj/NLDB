
public class Taylor {
	
	
	public static int factioal(int par1) {
		int ans = 1;
		for(int i =1; i<par1; i++) {
			ans *=i;
		}
		
		return ans;
	}
	/**
	 * 
	 * @param par1 - points set
	 * @param par2 - null at beginning
	 * @param i = 0 at beginning
	 * @return
	 */
	public static double[] generateFromPoints(double[] par1, double[] par2, int i) {
		if(par2 == null || par2.length == 0) {
			par2 = new double[par1.length];
			
			if(par1.length % 2 == 1) {
				par2[0] = par1[(par1.length / 2)];
			}
			if(par1.length % 2 == 0) {
				double sum =  par1[(par1.length / 2)] + par1[(par1.length / 2) - 1];
				par2[0] = sum/2;
			}
			
			
		
			i++;
			return generateFromPoints(par1, par2, i);
		}
		if(par1.length == 0 || par1.length == 1) {
			return par2;
		}
		
		
		int mod = (par1.length-1) % 2;
		double[] par3 = new double[par1.length-1];
		for(int j = 0; j < par3.length; j++) {
			par3[j] = par1[j+1]-par1[j];
		}
		if(mod == 1) {
			
			par2[i] = par3[(par3.length / 2)];
			i++;
		}
		else {
			double sum = par3[(par3.length / 2) - 1] + par3[(par3.length / 2)];//might have some errors on this line need to test this before i send it out
			par2[i] = sum/2;
			i++;
		}
		
		return generateFromPoints(par3, par2, i);
	}
	
	/**
	 * gets the distance squared for the line
	 * @param par1 data set1
	 * @param par2 data set2
	 * @return the distance squared
	 */
	public static double DistanceSquared(double[] par1, double[] par2) {
		double d = 0;
		for(int i =0; i < par1.length; i++) {
			d += Math.pow(par1[i] - par2[i], 2);
		}
		
		return d;
	}
	
	
	/**
	 * gets the distance squared for each point
	 * @param par1 data set1
	 * @param par2 data set2
	 * @return the distance squared
	 */
	public static double[] DistanceSquaredPoints(double[] par1, double[] par2) {
		double[] d = new double[par1.length];
		for(int i =0; i < par1.length; i++) {
			d[i]= Math.pow(par1[i] - par2[i], 2);
		}
		
		return d;
	}
	/**
	 * Generates Points from Taylor Series
	 * @param t Taylor Series
	 * @param par1 the orginal data
	 * @return points
	 */
	public static double[] dataSetFromTaylorSeries(TaylorSeries t) {
		double[] d = new double[t.getCoeff().length];
		
		for(int i = 0; i < t.getCoeff().length; i++) {
			d[i] = t.evaluteAt(i, (t.getCoeff().length-1) / 2.0);
		}
		
		return d;
	}
	
}
