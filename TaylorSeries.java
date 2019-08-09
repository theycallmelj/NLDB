
public class TaylorSeries {

	private double[] coeff;
	/**
	 * 
	 * @param par1 coefficients
	 */
	public TaylorSeries(double[] par1) {
		
			this.setCoeff(par1);
		
		
	}
	
	public double[] getCoeff() {
		return coeff;
	}
	public void setCoeff(double[] coeff) {
		this.coeff = coeff;
	}
	
	/**
	 * Evaluates the talyor series at x
	 * @param x the position
	 * @param a the center based on index
	 * @return the number
	 */
	public double evaluteAt(double x, double a) {
		double ans = 0.0;
		
		for(int i = 0; i < this.coeff.length; i++) {
			
			ans += Math.pow(x-a, i) * this.coeff[i]/ Taylor.factioal(i);
		}
		return ans;
	}

	
	
}
