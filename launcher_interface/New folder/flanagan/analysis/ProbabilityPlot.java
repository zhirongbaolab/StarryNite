/*
*   Class   ProbabilityPlot
*
*   USAGE:  Probability Plots
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   DATE:    29-30 September 2008, 1-5 October 2008
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web page:
*   http://www.ee.ucl.ac.uk/~mflanaga/java/ProbabilityPlot.html
*   http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*   Copyright (c) 2008 Michael Thomas Flanagan
*
*   PERMISSION TO COPY:
*
* Permission to use, copy and modify this software and its documentation for NON-COMMERCIAL purposes is granted, without fee,
* provided that an acknowledgement to the author, Dr Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies
* and associated documentation or publications.
*
* Redistributions of the source code of this source code, or parts of the source codes, must retain the above copyright notice, this list of conditions
* and the following disclaimer and requires written permission from the Michael Thomas Flanagan:
*
* Redistribution in binary form of all or parts of this class must reproduce the above copyright notice, this list of conditions and
* the following disclaimer in the documentation and/or other materials provided with the distribution and requires written permission from the Michael Thomas Flanagan:
*
* Dr Michael Thomas Flanagan makes no representations about the suitability or fitness of the software for any or for a particular purpose.
* Dr Michael Thomas Flanagan shall not be liable for any damages suffered as a result of using, modifying or distributing this software
* or its derivatives.
*
***************************************************************************************/

package flanagan.analysis;

import java.util.*;
import java.math.*;

import flanagan.math.*;
import flanagan.plot.PlotGraph;
import flanagan.analysis.*;
import flanagan.io.PrintToScreen;

public class ProbabilityPlot{

        // INSTANCE VARIABLES
        private double[] array = null;                      // array of data
        private Stat arrayAsStat = null;                    // array of data as Stat
        private double[] sortedData = null;                 // data sorted into ascending order

        private double mean = Double.NaN;                   // array mean
        private double standardDeviation = Double.NaN;      // array standard deviation
        private double minimum = Double.NaN;                // array minimum
        private double maximum = Double.NaN;                // array maximum
        private double range = Double.NaN;                  // array range

        private int numberOfDataPoints = 0;                 // number of data points

        private int gaussianNumberOfParameters = 2;         // number of Gaussian parameters
        private double[] gaussianOrderMedians = null;       // Gaussian order statistic medians
        private double[] gaussianParam = null;              // Gaussian parameters obtained by the minimization procedure
        private double[] gaussianParamErrors = null;        // estimates of the errors of the Gaussian parameters obtained by the minimization procedure
        private double gaussianSumOfSquares = Double.NaN;   // sum of squares at Gaussian minimum
        private double[] gaussianLine = null;               // Gaussian probability plot gradient and intercept
        private double[] gaussianLineErrors = null;         // estimated errors of the Gaussian probability plot gradient and intercept
        private double gaussianCorrCoeff = Double.NaN;      // Gaussian correlation coefficient of the probability plot
        private boolean gaussianDone = false;               // = true after Gaussian probability plot drawn

        private int weibullNumberOfParameters = 3;          // number of Weibull parameters
        private double[] weibullOrderMedians = null;        // Weibull order statistic medians
        private double[] weibullParam = null;               // Weibull parameters obtained by the minimization procedure
        private double[] weibullParamErrors = null;         // estimates of the errors of the Weibull parameters obtained by the minimization procedure
        private double weibullSumOfSquares = Double.NaN;    // sum of squares at Weibull minimum
        private double[] weibullLine = null;                // Weibull probability plot gradient and intercept
        private double[] weibullLineErrors = null;          // estimated errors of the Weibull probability plot gradient and intercept
        private double weibullCorrCoeff = Double.NaN;       // Weibull correlation coefficient of the probability plot
        private boolean weibullDone = false;                // = true after Weibull probability plot drawn

        private int exponentialNumberOfParameters = 2;      // number of Exponential parameters
        private double[] exponentialOrderMedians = null;    // Exponential order statistic medians
        private double[] exponentialParam = null;           // Exponential parameters obtained by the minimization procedure
        private double[] exponentialParamErrors = null;     // estimates of the errors of the Exponential parameters obtained by the minimization procedure
        private double exponentialSumOfSquares = Double.NaN;// sum of squares at Exponential minimum
        private double[] exponentialLine = null;            // Exponential probability plot gradient and intercept
        private double[] exponentialLineErrors = null;      // estimated errors of the Exponential probability plot gradient and intercept
        private double exponentialCorrCoeff = Double.NaN;   // Exponential correlation coefficient of the probability plot
        private boolean exponentialDone = false;            // = true after Exponential probability plot drawn

        private int rayleighNumberOfParameters = 2;         // number of Rayleigh parameters
        private double[] rayleighOrderMedians = null;       // Rayleigh order statistic medians
        private double[] rayleighParam = null;              // Rayleigh parameters obtained by the minimization procedure
        private double[] rayleighParamErrors = null;        // estimates of the errors of the Rayleigh parameters obtained by the minimization procedure
        private double rayleighSumOfSquares = Double.NaN;   // sum of squares at Rayleigh minimum
        private double[] rayleighLine = null;               // Rayleigh probability plot gradient and intercept
        private double[] rayleighLineErrors = null;         // estimated errors of the Rayleigh probability plot gradient and intercept
        private double rayleighCorrCoeff = Double.NaN;      // Rayleigh correlation coefficient of the probability plot
        private boolean rayleighDone = false;               // = true after Rayleigh probability plot drawn

        private int paretoNumberOfParameters = 2;           // number of Pareto parameters
        private double[] paretoOrderMedians = null;         // Pareto order statistic medians
        private double[] paretoParam = null;                // Pareto parameters obtained by the minimization procedure
        private double[] paretoParamErrors = null;          // estimates of the errors of the Pareto parameters obtained by the minimization procedure
        private double paretoSumOfSquares = Double.NaN;     // sum of squares at Pareto minimum
        private double[] paretoLine = null;                 // Pareto probability plot gradient and intercept
        private double[] paretoLineErrors = null;           // estimated errors of the Pareto probability plot gradient and intercept
        private double paretoCorrCoeff = Double.NaN;        // Pareto correlation coefficient of the probability plot
        private boolean paretoDone = false;                 // = true after Pareto probability plot drawn

        private boolean probPlotDone = false;               // = true after any probability plot drawn

        private double delta = 1e2;                        // step fraction in numerical differentiation

        private boolean nFactorOptionI = false;             // = true  variance, covariance and standard deviation denominator = n                                                                                // = false varaiance, covariance and standard deviation denominator = n-1
        private boolean nFactorReset = false;               // = true when instance method resetting the denominator is called



        // CONSTRUCTORS
        public ProbabilityPlot(double[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Double[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(float[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Float[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(long[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Long[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(int[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Integer[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(short[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Short[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(byte[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Byte[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(BigDecimal[] xx){
             this.arrayAsStat = new Stat(xx);
             this.initialize();
        }

        public ProbabilityPlot(BigInteger[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Object[] xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Vector<Object> xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(ArrayList<Object> xx){
             this.arrayAsStat = new Stat(xx);
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(ArrayMaths xx){
             this.arrayAsStat = xx.toStat();
             this.array = this.arrayAsStat.array();
             this.initialize();
        }

        public ProbabilityPlot(Stat xx){
             this.arrayAsStat = xx;
             this.array = this.arrayAsStat.array();
             this.initialize();
        }


        // INITIALIZATIONS
        private void initialize(){
            this.numberOfDataPoints = this.array.length;
            Stat sorted = arrayAsStat.sort();
            this.sortedData = sorted.array();
            this.mean = arrayAsStat.mean();
            this.minimum = arrayAsStat.minimum();
            this.maximum = arrayAsStat.maximum();
            this.range = this.maximum - this.minimum;
        }



        // GAUSSIAN PROBABILITY PLOT
        public void gaussianProbabilityPlot(){

            if(this.nFactorOptionI)arrayAsStat.setDenominatorToN();
            this.standardDeviation = arrayAsStat.standardDeviation();
            this.gaussianNumberOfParameters = 2;
            if(this.numberOfDataPoints<(this.gaussianNumberOfParameters+1))throw new IllegalArgumentException("There must be at least three data points - preferably considerably more");

            // Create instance of Minimization
            Minimization min = new Minimization();
            double meanest = this.mean;
            if(this.mean==0)meanest = this.standardDeviation/3.0;
            double[] start = {meanest, this.standardDeviation};
            double[] step = {0.3*meanest, 0.3*this.standardDeviation};
            double tolerance = 1e-4;

            // Add constraint; sigma>0
            min.addConstraint(1, -1, 0);

            // Create an instance of GaussProbPlotFunc
            GaussProbPlotFunc gppf = new GaussProbPlotFunc();
            gppf.setDataArray(sortedData);

            // Obtain best probability plot varying mu and sigma
            // by minimizing the sum of squares of the differences between the ordered data and the ordered statistic medians
            min.nelderMead(gppf, start, step, tolerance);

            // Get mu and sigma for best correlation coefficient
            this.gaussianParam = min.getParamValues();

            // Calculate Gaussian order statistic medians
            this.gaussianOrderMedians = Stat.gaussianOrderStatisticMedians(this.gaussianParam[0], this.gaussianParam[1], this.numberOfDataPoints);

            // Regression of the ordered data on the Gaussian order statistic medians
            Regression reg = new Regression(this.gaussianOrderMedians, this.sortedData);
            reg.linear();

            // Intercept and gradient of best fit straight line
            this.gaussianLine = reg.getBestEstimates();

            // Estimated erors of the intercept and gradient of best fit straight line
            this.gaussianLineErrors = reg.getBestEstimatesErrors();

            // Correlation coefficient
            this.gaussianCorrCoeff = reg.getSampleR();

            // Initialize data arrays for plotting
            double[][] data = PlotGraph.data(2,this.numberOfDataPoints);

            // Assign data to plotting arrays
            data[0] = this.gaussianOrderMedians;
            data[1] = this.sortedData;

            data[2] = this.gaussianOrderMedians;
            for(int i=0; i<this.numberOfDataPoints; i++){
                data[3][i] = this.gaussianLine[0] + this.gaussianLine[1]*this.gaussianOrderMedians[i];
            }

            // Estimates of the errors in the best fit mean and sd
            this.probPlotStats(0, gaussianParam);

            // Create instance of PlotGraph
            PlotGraph pg = new PlotGraph(data);
            int[] points = {4, 0};
            pg.setPoint(points);
            int[] lines = {0, 3};
            pg.setLine(lines);
            pg.setXaxisLegend("Gaussian Order Statistic Medians");
            pg.setYaxisLegend("Ordered Data Values");
            pg.setGraphTitle("Gaussian probability plot:   gradient = " + Fmath.truncate(this.gaussianLine[1], 4) + ", intercept = "  +  Fmath.truncate(this.gaussianLine[0], 4) + ",  R = " + Fmath.truncate(this.gaussianCorrCoeff, 4));
            pg.setGraphTitle2("  mu = " + Fmath.truncate(this.gaussianParam[0], 4) + ", sigma = "  +  Fmath.truncate(this.gaussianParam[1], 4));

            // Plot
            pg.plot();

            this.gaussianDone = true;
        }

        public void normalProbabilityPlot(){
            this.gaussianProbabilityPlot();
        }

        // Return Gaussian mu
        public double gaussianMu(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParam[0];
        }

        // Return Gaussian mu error
        public double gaussianMuError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParamErrors[0];
        }

        // Return Gaussian sigma
        public double gaussianSigma(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParam[1];
        }

        // Return Gaussian sigma error
        public double gaussianSigmaError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParamErrors[1];
        }

        // Return the Gaussian gradient
        public double gaussianGradient(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLine[1];
        }

        // Return the error of the Gaussian gradient
        public double gaussianGradientError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLineErrors[1];
        }

        // Return the Gaussian intercept
        public double gaussianIntercept(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLine[0];
        }

        // Return the error of the Gaussian intercept
        public double gaussianInterceptError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLineErrors[0];
        }

        // Return the Gaussian correlation coefficient
        public double gaussianCorrelationCoefficient(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianCorrCoeff;
        }

        // Return the sum of squares at the Gaussian minimum
        public double gaussianSumOfSquares(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianSumOfSquares;
        }

        // Return Gaussian order statistic medians
        public double[] gaussianOrderStatisticMedians(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianOrderMedians;
        }


        public double normalMu(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParam[0];
        }

        // Return Gaussian mu error
        public double normalMuError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParamErrors[0];
        }

        // Return Gaussian sigma
        public double normalSigma(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParam[1];
        }

        // Return Gaussian sigma error
        public double normalSigmaError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianParamErrors[1];
        }

        // Return the Gaussian gradient
        public double normalGradient(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLine[1];
        }

        // Return the error of the Gaussian gradient
        public double normalGradientError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLineErrors[1];
        }

        // Return the Gaussian intercept
        public double normalIntercept(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLine[0];
        }

        // Return the error of the Gaussian intercept
        public double normalInterceptError(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianLineErrors[0];
        }

        // Return the Gaussian correlation coefficient
        public double normalCorrelationCoefficient(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianCorrCoeff;
        }

        // Return the sum of squares at the Gaussian minimum
        public double normalSumOfSquares(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianSumOfSquares;
        }

        // Return Gaussian order statistic medians
        public double[] normalOrderStatisticMedians(){
            if(!this.gaussianDone)throw new IllegalArgumentException("Gaussian Probability Plot method has not been called");
            return this.gaussianOrderMedians;
        }





        // WEIBULL PROBABILITY PLOT
        public void weibullProbabilityPlot(){

            if(this.nFactorOptionI)arrayAsStat.setDenominatorToN();
            this.standardDeviation = arrayAsStat.standardDeviation();
            this.weibullNumberOfParameters = 3;
            if(this.numberOfDataPoints<4)throw new IllegalArgumentException("There must be at least four data points - preferably considerably more");

            // Create instance of Minimization
            Minimization min = new Minimization();
            double muest = minimum;
            if(muest==0.0)muest = this.standardDeviation/3.0;
            double sigmaest = this.standardDeviation;
            double gammaest = 4.0;
            double[] start = {muest, sigmaest, gammaest};
            double[] step = {0.3*muest, 0.3*sigmaest, 0.3*gammaest};
            double tolerance = 1e-4;

             // Add constraint; sigma>0, gamma>0
            min.addConstraint(1, -1, 0);
            min.addConstraint(2, -1, 0);

            // Create an instance of WeibullProbPlotFunc
            WeibullProbPlotFunc wppf = new WeibullProbPlotFunc();
            wppf.setDataArray(this.sortedData);

            // Obtain best probability plot varying mu and sigma
            // by minimizing the sum of squares of the differences between the ordered data and the ordered statistic medians
            min.nelderMead(wppf, start, step, tolerance);

            // Get mu, sigma and gamma for best correlation coefficient
            this.weibullParam = min.getParamValues();

            // Calculate Weibull order statistic medians
            this.weibullOrderMedians = Stat.weibullOrderStatisticMedians(this.weibullParam[0], this.weibullParam[1], this.weibullParam[2], this.numberOfDataPoints);

            // Regression of the ordered data on the Weibull order statistic medians
            Regression reg = new Regression(this.weibullOrderMedians, this.sortedData);
            reg.linear();

            // Intercept and gradient of best fit straight line
            this.weibullLine = reg.getBestEstimates();

            // Estimated erors of the intercept and gradient of best fit straight line
            this.weibullLineErrors = reg.getBestEstimatesErrors();

            // Correlation coefficient
            this.weibullCorrCoeff = reg.getSampleR();

            // Initialize data arrays for plotting
            double[][] data = PlotGraph.data(2,this.numberOfDataPoints);

            // Assign data to plotting arrays
            data[0] = this.weibullOrderMedians;
            data[1] = this.sortedData;

            data[2] = weibullOrderMedians;
            for(int i=0; i<this.numberOfDataPoints; i++){
                data[3][i] = this.weibullLine[0] + this.weibullLine[1]*weibullOrderMedians[i];
            }

            // Estimates of the errors in the best fit mu, sigma and gamma
            this.probPlotStats(1, this.weibullParam);

            // Create instance of PlotGraph
            PlotGraph pg = new PlotGraph(data);
            int[] points = {4, 0};
            pg.setPoint(points);
            int[] lines = {0, 3};
            pg.setLine(lines);
            pg.setXaxisLegend("Weibull Order Statistic Medians");
            pg.setYaxisLegend("Ordered Data Values");
            pg.setGraphTitle("Weibull probability plot:   gradient = " + Fmath.truncate(this.weibullLine[1], 4) + ", intercept = "  +  Fmath.truncate(this.weibullLine[0], 4) + ",  R = " + Fmath.truncate(this.weibullCorrCoeff, 4));
            pg.setGraphTitle2("  mu = " + Fmath.truncate(this.weibullParam[0], 4) + ", sigma = "  +  Fmath.truncate(this.weibullParam[1], 4) + ", gamma = "  +  Fmath.truncate(this.weibullParam[2], 4));

            // Plot
            pg.plot();

            this.weibullDone = true;
            this.probPlotDone = true;
        }

        // Return Weibull mu
        public double weibullMu(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParam[0];
        }

        // Return Weibull mu error
        public double weibullMuError(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParamErrors[0];
        }

        // Return Weibull sigma
        public double weibullSigma(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParam[1];
        }

        // Return Weibull sigma error
        public double weibullSigmaError(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParamErrors[1];
        }

        // Return Weibull gamma
        public double weibullGamma(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParam[2];
        }

        // Return Weibull gamma error
        public double weibullGammaError(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullParamErrors[2];
        }

        // Return Weibull order statistic medians
        public double[] weibullOrderStatisticMedians(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullOrderMedians;
        }

        // Return the Weibull gradient
        public double weibullGradient(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullLine[1];
        }

        // Return the error of the Weibull gradient
        public double weibullGradientError(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullLineErrors[1];
        }

        // Return the Weibull intercept
        public double weibullIntercept(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullLine[0];
        }

        // Return the error of the Weibull intercept
        public double weibullInterceptError(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullLineErrors[0];
        }

        // Return the Weibull correlation coefficient
        public double weibullCorrelationCoefficient(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullCorrCoeff;
        }

        // Return the sum of squares at the Weibull minimum
        public double weibullSumOfSquares(){
            if(!this.weibullDone)throw new IllegalArgumentException("Weibull Probability Plot method has not been called");
            return this.weibullSumOfSquares;
        }

        // EXPONENTIAL PROBABILITY PLOT
        public void exponentialProbabilityPlot(){

            if(this.nFactorOptionI)arrayAsStat.setDenominatorToN();
            this.standardDeviation = arrayAsStat.standardDeviation();
            this.exponentialNumberOfParameters = 2;
            if(this.numberOfDataPoints<3)throw new IllegalArgumentException("There must be at least three data points - preferably considerably more");

            // Create instance of Minimization
            Minimization min = new Minimization();
            double muest = minimum;
            if(muest==0.0)muest = this.standardDeviation/3.0;
            double sigmaest = this.standardDeviation;
            double[] start = {muest, sigmaest};
            double[] step = {0.3*muest, 0.3*sigmaest};
            double tolerance = 1e-4;

             // Add constraint; sigma>0
            min.addConstraint(1, -1, 0);

            // Create an instance of ExponentialProbPlotFunc
            ExponentialProbPlotFunc eppf = new ExponentialProbPlotFunc();
            eppf.setDataArray(this.sortedData);

            // Obtain best probability plot varying mu and sigma
            // by minimizing the sum of squares of the differences between the ordered data and the ordered statistic medians
            min.nelderMead(eppf, start, step, tolerance);

            // Get mu, sigma and gamma values
            this.exponentialParam = min.getParamValues();

            // Calculate Exponential order statistic medians (Weibull with gamma = 1)
            this.exponentialOrderMedians = Stat.weibullOrderStatisticMedians(this.exponentialParam[0], this.exponentialParam[1], 1.0, this.numberOfDataPoints);

            // Regression of the ordered data on the Exponential order statistic medians
            Regression reg = new Regression(this.exponentialOrderMedians, this.sortedData);
            reg.linear();

            // Intercept and gradient of best fit straight line
            this.exponentialLine = reg.getBestEstimates();

            // Estimated erors of the intercept and gradient of best fit straight line
            this.exponentialLineErrors = reg.getBestEstimatesErrors();

            // Correlation coefficient
            this.exponentialCorrCoeff = reg.getSampleR();

            // Initialize data arrays for plotting
            double[][] data = PlotGraph.data(2,this.numberOfDataPoints);

            // Assign data to plotting arrays
            data[0] = this.exponentialOrderMedians;
            data[1] = this.sortedData;

            data[2] = exponentialOrderMedians;
            for(int i=0; i<this.numberOfDataPoints; i++){
                data[3][i] = this.exponentialLine[0] + this.exponentialLine[1]*exponentialOrderMedians[i];
            }

            // Estimates of the errors in the best fit mu, sigma and gamma
            this.probPlotStats(2, this.exponentialParam);

            // Create instance of PlotGraph
            PlotGraph pg = new PlotGraph(data);
            int[] points = {4, 0};
            pg.setPoint(points);
            int[] lines = {0, 3};
            pg.setLine(lines);
            pg.setXaxisLegend("Exponential Order Statistic Medians");
            pg.setYaxisLegend("Ordered Data Values");
            pg.setGraphTitle("Exponential probability plot:   gradient = " + Fmath.truncate(this.exponentialLine[1], 4) + ", intercept = "  +  Fmath.truncate(this.exponentialLine[0], 4) + ",  R = " + Fmath.truncate(this.exponentialCorrCoeff, 4));
            pg.setGraphTitle2("  mu = " + Fmath.truncate(this.exponentialParam[0], 4) + ", sigma = "  +  Fmath.truncate(this.exponentialParam[1], 4));

            // Plot
            pg.plot();

            this.exponentialDone = true;
            this.probPlotDone = true;
        }

        // Return Exponential mu
        public double exponentialMu(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialParam[0];
        }

        // Return Exponential mu error
        public double exponentialMuError(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialParamErrors[0];
        }

        // Return Exponential sigma
        public double exponentialSigma(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialParam[1];
        }

        // Return Exponential sigma error
        public double exponentialSigmaError(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialParamErrors[1];
        }


        // Return Exponential order statistic medians
        public double[] exponentialOrderStatisticMedians(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialOrderMedians;
        }

        // Return the Exponential gradient
        public double exponentialGradient(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialLine[1];
        }

        // Return the error of the Exponential gradient
        public double exponentialGradientError(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialLineErrors[1];
        }

        // Return the Exponential intercept
        public double exponentialIntercept(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialLine[0];
        }

        // Return the error of the Exponential intercept
        public double exponentialInterceptError(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialLineErrors[0];
        }

        // Return the Exponential correlation coefficient
        public double exponentialCorrelationCoefficient(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialCorrCoeff;
        }

        // Return the sum of squares at the Exponential minimum
        public double exponentialSumOfSquares(){
            if(!this.exponentialDone)throw new IllegalArgumentException("Exponential Probability Plot method has not been called");
            return this.exponentialSumOfSquares;
        }


        // RAYLEIGH PROBABILITY PLOT
        public void rayleighProbabilityPlot(){

            if(this.nFactorOptionI)arrayAsStat.setDenominatorToN();
            this.standardDeviation = arrayAsStat.standardDeviation();
            this.rayleighNumberOfParameters = 1;
            if(this.numberOfDataPoints<2)throw new IllegalArgumentException("There must be at least two data points - preferably considerably more");

            // Create instance of Minimization
            Minimization min = new Minimization();
             double sigmaest = this.standardDeviation;
            double[] start = {sigmaest};
            double[] step = {0.3*sigmaest};
            double tolerance = 1e-4;

             // Add constraint; beta>0
            min.addConstraint(0, -1, 0);

            // Create an instance of RayleighProbPlotFunc
            RayleighProbPlotFunc rppf = new RayleighProbPlotFunc();
            rppf.setDataArray(this.sortedData);

            // Obtain best probability plot varying mu and sigma
            // by minimizing the sum of squares of the differences between the ordered data and the ordered statistic medians
            min.nelderMead(rppf, start, step, tolerance);

            // Get mu, sigma and gamma values
            this.rayleighParam = min.getParamValues();

            // Calculate Rayleigh order statistic medians (Weibull with mu = 0, sigma = sqrt(2).beta, gamma = 2)
            this.rayleighOrderMedians = Stat.weibullOrderStatisticMedians(0.0, this.rayleighParam[0]*Math.sqrt(2.0), 2.0, this.numberOfDataPoints);

            // Regression of the ordered data on the Rayleigh order statistic medians
            Regression reg = new Regression(this.rayleighOrderMedians, this.sortedData);
            reg.linear();

            // Intercept and gradient of best fit straight line
            this.rayleighLine = reg.getBestEstimates();

            // Estimated erors of the intercept and gradient of best fit straight line
            this.rayleighLineErrors = reg.getBestEstimatesErrors();

            // Correlation coefficient
            this.rayleighCorrCoeff = reg.getSampleR();

            // Initialize data arrays for plotting
            double[][] data = PlotGraph.data(2,this.numberOfDataPoints);

            // Assign data to plotting arrays
            data[0] = this.rayleighOrderMedians;
            data[1] = this.sortedData;

            data[2] = rayleighOrderMedians;
            for(int i=0; i<this.numberOfDataPoints; i++){
                data[3][i] = this.rayleighLine[0] + this.rayleighLine[1]*rayleighOrderMedians[i];
            }

            // Estimates of the errors in the best fit mu, sigma and gamma
            this.probPlotStats(3, this.rayleighParam);

            // Create instance of PlotGraph
            PlotGraph pg = new PlotGraph(data);
            int[] points = {4, 0};
            pg.setPoint(points);
            int[] lines = {0, 3};
            pg.setLine(lines);
            pg.setXaxisLegend("Rayleigh Order Statistic Medians");
            pg.setYaxisLegend("Ordered Data Values");
            pg.setGraphTitle("Rayleigh probability plot:   gradient = " + Fmath.truncate(this.rayleighLine[1], 4) + ", intercept = "  +  Fmath.truncate(this.rayleighLine[0], 4) + ",  R = " + Fmath.truncate(this.rayleighCorrCoeff, 4));
            pg.setGraphTitle2("  beta = " + Fmath.truncate(this.rayleighParam[0], 4));

            // Plot
            pg.plot();

            this.rayleighDone = true;
            this.probPlotDone = true;
        }

        // Return Rayleigh beta
        public double rayleighBeta(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighParam[0];
        }

        // Return Rayleigh beta error
        public double rayleighBetaError(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighParamErrors[0];
        }

        // Return Rayleigh order statistic medians
        public double[] rayleighOrderStatisticMedians(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighOrderMedians;
        }

        // Return the Rayleigh gradient
        public double rayleighGradient(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighLine[1];
        }

        // Return the error of the Rayleigh gradient
        public double rayleighGradientError(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighLineErrors[1];
        }

        // Return the Rayleigh intercept
        public double rayleighIntercept(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighLine[0];
        }

        // Return the error of the Rayleigh intercept
        public double rayleighInterceptError(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighLineErrors[0];
        }

        // Return the Rayleigh correlation coefficient
        public double rayleighCorrelationCoefficient(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighCorrCoeff;
        }

        // Return the sum of squares at the Rayleigh minimum
        public double rayleighSumOfSquares(){
            if(!this.rayleighDone)throw new IllegalArgumentException("Rayleigh Probability Plot method has not been called");
            return this.rayleighSumOfSquares;
        }

        // PARETO PROBABILITY PLOT
        public void paretoProbabilityPlot(){

            if(this.nFactorOptionI)arrayAsStat.setDenominatorToN();
            this.standardDeviation = arrayAsStat.standardDeviation();
            this.paretoNumberOfParameters = 2;
            if(this.numberOfDataPoints<3)throw new IllegalArgumentException("There must be at least two data points - preferably considerably more");

            // Create instance of Minimization
            Minimization min = new Minimization();
            double betaest = this.minimum;
            double alphaest = this.mean/(this.mean - betaest);
            double[] start = {alphaest, betaest};
            double[] step = {0.3*alphaest, 0.3*betaest};
            double tolerance = 1e-4;

            // Create an instance of ParetoProbPlotFunc
            ParetoProbPlotFunc pppf = new ParetoProbPlotFunc();
            pppf.setDataArray(this.sortedData);

            // Obtain best probability plot varying mu and sigma
            // by minimizing the sum of squares of the differences between the ordered data and the ordered statistic medians
            min.nelderMead(pppf, start, step, tolerance);

            // Get mu, sigma and gamma values
            this.paretoParam = min.getParamValues();

            // Calculate Pareto order statistic medians
            this.paretoOrderMedians = Stat.paretoOrderStatisticMedians(this.paretoParam[0], this.paretoParam[1], this.numberOfDataPoints);

            // Regression of the ordered data on the Pareto order statistic medians
            Regression reg = new Regression(this.paretoOrderMedians, this.sortedData);
            reg.linear();

            // Intercept and gradient of best fit straight line
            this.paretoLine = reg.getBestEstimates();

            // Estimated erors of the intercept and gradient of best fit straight line
            this.paretoLineErrors = reg.getBestEstimatesErrors();

            // Correlation coefficient
            this.paretoCorrCoeff = reg.getSampleR();

            // Initialize data arrays for plotting
            double[][] data = PlotGraph.data(2,this.numberOfDataPoints);

            // Assign data to plotting arrays
            data[0] = this.paretoOrderMedians;
            data[1] = this.sortedData;

            data[2] = paretoOrderMedians;
            for(int i=0; i<this.numberOfDataPoints; i++){
                data[3][i] = this.paretoLine[0] + this.paretoLine[1]*paretoOrderMedians[i];
            }

            // Estimates of the errors in the best fit mu, sigma and gamma
            this.probPlotStats(4, this.paretoParam);

            // Create instance of PlotGraph
            PlotGraph pg = new PlotGraph(data);
            int[] points = {4, 0};
            pg.setPoint(points);
            int[] lines = {0, 3};
            pg.setLine(lines);
            pg.setXaxisLegend("Pareto Order Statistic Medians");
            pg.setYaxisLegend("Ordered Data Values");
            pg.setGraphTitle("Pareto probability plot:   gradient = " + Fmath.truncate(this.paretoLine[1], 4) + ", intercept = "  +  Fmath.truncate(this.paretoLine[0], 4) + ",  R = " + Fmath.truncate(this.paretoCorrCoeff, 4));
            pg.setGraphTitle2("  alpha = " + Fmath.truncate(this.paretoParam[0], 4) + ", beta = "  +  Fmath.truncate(this.paretoParam[1], 4));

            // Plot
            pg.plot();

            this.paretoDone = true;
            this.probPlotDone = true;
        }

        // Return Pareto alpha
        public double paretoAlpha(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoParam[0];
        }

        // Return Pareto alpha error
        public double paretoAlphaError(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoParamErrors[0];
        }

        // Return Pareto beta
        public double paretoBeta(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoParam[1];
        }

        // Return Pareto beta error
        public double paretoBetaError(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoParamErrors[1];
        }

        // Return Pareto order statistic medians
        public double[] paretoOrderStatisticMedians(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoOrderMedians;
        }

        // Return the Pareto gradient
        public double paretoGradient(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoLine[1];
        }

        // Return the error of the Pareto gradient
        public double paretoGradientError(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoLineErrors[1];
        }

        // Return the Pareto intercept
        public double paretoIntercept(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoLine[0];
        }

        // Return the error of the Pareto intercept
        public double paretoInterceptError(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoLineErrors[0];
        }

        // Return the Pareto correlation coefficient
        public double paretoCorrelationCoefficient(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoCorrCoeff;
        }

        // Return the sum of squares at the Pareto minimum
        public double paretoSumOfSquares(){
            if(!this.paretoDone)throw new IllegalArgumentException("Pareto Probability Plot method has not been called");
            return this.paretoSumOfSquares;
        }




        // COMMON METHODS

        // Return the ordered data
        public double[] orderedData(){
            return this.sortedData;
        }

        // Return the number of data points
        public int numberOfDataPoints(){
            return this.numberOfDataPoints;
        }

        // Return the data mean
        public double mean(){
            return this.mean;
        }

        // Return the data standard deviation
        public double standardDeviation(){
            if(!this.probPlotDone)throw new IllegalArgumentException("no probability plot method has been called");
            return this.standardDeviation;
        }

        // Return the data minimum
        public double minimum(){
            return this.minimum;
        }

        // Return the data maximum
        public double maximum(){
             return this.maximum;
        }

        // Return the numerical differentiation step, delta
        public double delta(){
            return this.delta;
        }

        // Reset the numerical differentiation step, delta
        public void resetDelta(double delta){
            this.delta = delta;
        }

        // Set standard deviation denominator to n
        public void setDenominatorToN(){
            this.nFactorOptionI = true;
            this.nFactorReset = true;
        }

        // Set standard deviation denominator to n-1
        public void setDenominatorToNminusOne(){
            this.nFactorOptionI = false;
            this.nFactorReset = true;
        }



        // STATISTICS
        // Estimate of the errors of the minimization parameters
        private void probPlotStats(int flag, double[] minParam){

            double range = this.maximum - this.minimum;
            double[][] secondDeriv = null;
            double[] paramTemp = null;
            int nMatrix = 0;
            GaussProbPlotFunc gppf = null;
            WeibullProbPlotFunc wppf = null;
            ExponentialProbPlotFunc eppf = null;
            RayleighProbPlotFunc rppf = null;
            ParetoProbPlotFunc pppf = null;
            Object func = null;
            switch(flag){
                case 0: //Gaussian Probability Plot
                        this.gaussianParamErrors = new double[this.gaussianNumberOfParameters];
                        gppf = new GaussProbPlotFunc();
                        secondDeriv = new double[this.gaussianNumberOfParameters][this.gaussianNumberOfParameters];
                        paramTemp = new double[this.gaussianNumberOfParameters];
                        gppf.setDataArray(this.array);
                        this.gaussianSumOfSquares = gppf.function(minParam);
                        func = (Object)gppf;
                        nMatrix = this.gaussianNumberOfParameters;
                        break;
                case 1: //Weibull Probability Plot
                        this.weibullParamErrors = new double[this.weibullNumberOfParameters];
                        wppf = new WeibullProbPlotFunc();
                        secondDeriv = new double[this.weibullNumberOfParameters][this.weibullNumberOfParameters];
                        paramTemp = new double[this.weibullNumberOfParameters];
                        wppf.setDataArray(this.array);
                        this.weibullSumOfSquares = wppf.function(minParam);
                        func = (Object)wppf;
                        nMatrix = this.weibullNumberOfParameters;
                        break;
                case 2: //Exponential Probability Plot
                        this.exponentialParamErrors = new double[this.exponentialNumberOfParameters];
                        eppf = new ExponentialProbPlotFunc();
                        secondDeriv = new double[this.exponentialNumberOfParameters][this.exponentialNumberOfParameters];
                        paramTemp = new double[this.exponentialNumberOfParameters];
                        eppf.setDataArray(this.array);
                        this.exponentialSumOfSquares = eppf.function(minParam);
                        func = (Object)eppf;
                        nMatrix = this.exponentialNumberOfParameters;
                        break;
                case 3: //Rayleigh Probability Plot
                        this.rayleighParamErrors = new double[rayleighNumberOfParameters];
                        rppf = new RayleighProbPlotFunc();
                        secondDeriv = new double[this.rayleighNumberOfParameters][this.rayleighNumberOfParameters];
                        paramTemp = new double[this.rayleighNumberOfParameters];
                        rppf.setDataArray(this.array);
                        this.rayleighSumOfSquares = rppf.function(minParam);
                        func = (Object)rppf;
                        nMatrix = this.rayleighNumberOfParameters;
                        break;
                case 4: //Pareto Probability Plot
                        this.paretoParamErrors = new double[paretoNumberOfParameters];
                        pppf = new ParetoProbPlotFunc();
                        secondDeriv = new double[this.paretoNumberOfParameters][this.paretoNumberOfParameters];
                        paramTemp = new double[this.paretoNumberOfParameters];
                        pppf.setDataArray(this.array);
                        this.paretoSumOfSquares = pppf.function(minParam);
                        func = (Object)pppf;
                        nMatrix = this.paretoNumberOfParameters;
                        break;
            }


            double pari = 0.0;
            double parj = 0.0;
            double[] parh = new double[2];
            for(int i=0; i<nMatrix; i++){
                    for(int j=0; j<nMatrix; j++){
                            if(i==j){
                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 + this.delta);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }
                                    double term1 = this.statFunction(flag, func, paramTemp);

                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i];
                                    pari = minParam[i];
                                    double term2 = this.statFunction(flag, func, paramTemp);

                                    double term3 = term2;

                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 - this.delta);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = -this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }
                                    double term4 = this.statFunction(flag, func, paramTemp);

                                    secondDeriv[i][j] = (term1 - term2 - term3 + term4)/(pari*pari*this.delta*this.delta);
                            }
                            else{
                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 + this.delta/2.0);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }

                                    paramTemp[j] = minParam[j]*(1.0 + this.delta/2.0);
                                    parj = minParam[j];
                                    if(minParam[j]==0.0){
                                        minParam[j] = this.range*this.delta/70.0;
                                        parj = this.range/35.0;
                                    }
                                    double term1 = this.statFunction(flag, func, paramTemp);

                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 - this.delta/2.0);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = -this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }

                                    paramTemp[j] = minParam[j]*(1.0 + this.delta/2.0);
                                    parj = minParam[j];
                                    if(minParam[j]==0.0){
                                        minParam[j] = this.range*this.delta/70.0;
                                        parj = this.range/35.0;
                                    }
                                    double term2 = this.statFunction(flag, func, paramTemp);

                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 + this.delta/2.0);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }

                                    paramTemp[j] = minParam[j]*(1.0 - this.delta/2.0);
                                    parj = minParam[j];
                                    if(minParam[j]==0.0){
                                        minParam[j] = -this.range*this.delta/70.0;
                                        parj = this.range/35.0;
                                    }
                                    double term3 = this.statFunction(flag, func, paramTemp);

                                    paramTemp = minParam.clone();
                                    paramTemp[i] = minParam[i]*(1.0 - this.delta/2.0);
                                    pari = minParam[i];
                                    if(minParam[i]==0.0){
                                        minParam[i] = -this.range*this.delta/70.0;
                                        pari = this.range/35.0;
                                    }

                                    paramTemp[j] = minParam[j]*(1.0 - this.delta/2.0);
                                    parj = minParam[j];
                                    if(minParam[j]==0.0){
                                        minParam[j] = -this.range*this.delta/70.0;
                                        parj = this.range/35.0;
                                    }
                                    double term4 = this.statFunction(flag, func, paramTemp);

                                    secondDeriv[i][j] = (term1 - term2 - term3 + term4)/(pari*parj*this.delta*this.delta);
                            }
                    }
            }
            Matrix mat = new Matrix(secondDeriv);
            mat = mat.inverse();
            switch(flag){
                case 0: mat = mat.times(this.gaussianSumOfSquares/(this.numberOfDataPoints - this.gaussianNumberOfParameters));
                        for(int i=0; i<this.gaussianNumberOfParameters; i++){
                            this.gaussianParamErrors[i] = Math.sqrt(mat.getElement(i,i));
                            if(Double.isNaN(this.gaussianParamErrors[i]))this.gaussianParamErrors[i] = Math.abs(this.gaussianSumOfSquares/(secondDeriv[i][i]*(this.numberOfDataPoints - this.gaussianNumberOfParameters)));
                        }
                        break;
                case 1: mat = mat.times(this.weibullSumOfSquares/(this.numberOfDataPoints - this.weibullNumberOfParameters));
                        for(int i=0; i<this.weibullNumberOfParameters; i++){
                            this.weibullParamErrors[i] = Math.sqrt(mat.getElement(i,i));
                            if(Double.isNaN(this.weibullParamErrors[i]))this.weibullParamErrors[i] = Math.abs(this.weibullSumOfSquares/(secondDeriv[i][i]*(this.numberOfDataPoints - this.weibullNumberOfParameters)));
                        }
                        break;
                case 2: mat = mat.times(this.exponentialSumOfSquares/(this.numberOfDataPoints - this.exponentialNumberOfParameters));
                        for(int i=0; i<this.exponentialNumberOfParameters; i++){
                            this.exponentialParamErrors[i] = Math.sqrt(mat.getElement(i,i));
                            if(Double.isNaN(this.exponentialParamErrors[i]))this.exponentialParamErrors[i] = Math.abs(this.exponentialSumOfSquares/(secondDeriv[i][i]*(this.numberOfDataPoints - this.exponentialNumberOfParameters)));
                        }
                        break;
                case 3: mat = mat.times(this.rayleighSumOfSquares/(this.numberOfDataPoints - this.rayleighNumberOfParameters));
                        for(int i=0; i<this.rayleighNumberOfParameters; i++){
                            this.rayleighParamErrors[i] = Math.sqrt(mat.getElement(i,i));
                            if(Double.isNaN(this.rayleighParamErrors[i]))this.rayleighParamErrors[i] = Math.abs(this.rayleighSumOfSquares/(secondDeriv[i][i]*(this.numberOfDataPoints - this.rayleighNumberOfParameters)));
                        }
                        break;
                case 4: mat = mat.times(this.paretoSumOfSquares/(this.numberOfDataPoints - this.paretoNumberOfParameters));
                        for(int i=0; i<this.paretoNumberOfParameters; i++){
                            this.paretoParamErrors[i] = Math.sqrt(mat.getElement(i,i));
                            if(Double.isNaN(this.paretoParamErrors[i]))this.paretoParamErrors[i] = Math.abs(this.paretoSumOfSquares/(secondDeriv[i][i]*(this.numberOfDataPoints - this.paretoNumberOfParameters)));
                        }
                        break;
            }
        }

        // Sum of squares function
        private double statFunction(int flag, Object func, double[] param){
            double term = 0.0;
            switch(flag){
                case 0: // Gaussian
                        term = ((GaussProbPlotFunc)func).function(param);
                        break;
                case 1: // Weibull
                        term = ((WeibullProbPlotFunc)func).function(param);
                        break;
                case 2: // Exponential
                        term = ((ExponentialProbPlotFunc)func).function(param);
                        break;
                case 3: // Rayleigh
                        term = ((RayleighProbPlotFunc)func).function(param);
                        break;
                case 4: // Pareto
                        term = ((ParetoProbPlotFunc)func).function(param);
                        break;
            }
            return term;
        }


}



// PROBABILITY PLOT FUNCTIONS
// Gaussian Probabilty plot function
class GaussProbPlotFunc implements MinimizationFunction{

    private double[] dataArray = null;
    private double[] medians = null;

    public double function(double[] x){

        // Calculate Gaussian order statistic medians
        medians = Stat.gaussianOrderStatisticMedians(x[0], x[1], this.dataArray.length);

        // Calculate sum of squares of medians - data
        double sum = 0.0;
        for(int i=0; i<this.dataArray.length; i++)sum += Fmath.square(medians[i] - dataArray[i]);

        return sum;

    }

    public void setDataArray(double[] array){
        this.dataArray = array;
    }
}


// Weibull Probabilty plot function
class WeibullProbPlotFunc implements MinimizationFunction{

    private double[] dataArray = null;
    private double[] medians = null;

    public double function(double[] x){

        // Calculate Weibull order statistic medians
        medians = Stat.weibullOrderStatisticMedians(x[0], x[1], x[2], this.dataArray.length);

        // Calculate sum of squares of medians - data
        double sum = 0.0;
        for(int i=0; i<this.dataArray.length; i++)sum += Fmath.square(medians[i] - dataArray[i]);

        return sum;

    }

    public void setDataArray(double[] array){
        this.dataArray = array;
    }
}

// Exponential Probabilty plot function
class ExponentialProbPlotFunc implements MinimizationFunction{

    private double[] dataArray = null;
    private double[] medians = null;

    public double function(double[] x){

        // Calculate Weibull order statistic medians
         medians = Stat.weibullOrderStatisticMedians(x[0], x[1], 1.0, this.dataArray.length);

        // Calculate sum of squares of medians - data
        double sum = 0.0;
        for(int i=0; i<this.dataArray.length; i++)sum += Fmath.square(medians[i] - dataArray[i]);

        return sum;
    }

    public void setDataArray(double[] array){
        this.dataArray = array;
    }
}

// Rayleigh Probabilty plot function
class RayleighProbPlotFunc implements MinimizationFunction{

    private double[] dataArray = null;
    private double[] medians = null;
    private double sqrt2 = Math.sqrt(2.0);

    public double function(double[] x){

        // Calculate Rayleigh order statistic medians (Weibull with mu = 0, sigma = sqrt(2).beta, gamma = 2)
         medians = Stat.weibullOrderStatisticMedians(0.0, sqrt2*x[0], 2.0, this.dataArray.length);

        // Calculate sum of squares of medians - data
        double sum = 0.0;
        for(int i=0; i<this.dataArray.length; i++)sum += Fmath.square(medians[i] - dataArray[i]);

        return sum;
    }

    public void setDataArray(double[] array){
        this.dataArray = array;
    }
}


// Pareto Probabilty plot function
class ParetoProbPlotFunc implements MinimizationFunction{

    private double[] dataArray = null;
    private double[] medians = null;

    public double function(double[] x){

        // Calculate Pareto order statistic medians
         medians = Stat.paretoOrderStatisticMedians(x[0], x[1], this.dataArray.length);

        // Calculate sum of squares of medians - data
        double sum = 0.0;
        for(int i=0; i<this.dataArray.length; i++)sum += Fmath.square(medians[i] - dataArray[i]);

        return sum;
    }

    public void setDataArray(double[] array){
        this.dataArray = array;
    }
}




