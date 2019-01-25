/*      Class ClosedLoop
*
*       This class supports the creation of a path of Black Boxes
*       i.e. of instances of BlackBox and of any of its subclasses,
*       e.g. PropIntDeriv, FirstOrder, and the methods to combine
*       these into both a single instance of BlackBox and a Vector
*       of analogue segments, digital segments and converters,
*       with a feedback path from the last box on the forward path to the first box on the forward path
*
*       Author:  Michael Thomas Flanagan.
*
*       Created: August 2002
*	    Updated: 14 May 2005, 6 April 2008, 5 July 2008
*
*       DOCUMENTATION:
*       See Michael T Flanagan's JAVA library on-line web page:
*       http://www.ee.ucl.ac.uk/~mflanaga/java/OpenLoop.html
*       http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*
*   Copyright (c) 2002 - 2008   Michael Thomas Flanagan
*
*   PERMISSION TO COPY:
*   Permission to use, copy and modify this software and its documentation for
*   NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*   to the author, Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*   Dr Michael Thomas Flanagan makes no representations about the suitability
*   or fitness of the software for any or for a particular purpose.
*   Michael Thomas Flanagan shall not be liable for any damages suffered
*   as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/

package flanagan.control;

import java.util.Vector;
import java.util.ArrayList;
import flanagan.complex.Complex;
import flanagan.complex.ComplexPoly;
import flanagan.control.OpenLoop;

public class ClosedLoop extends BlackBox{
    private OpenLoop forwardPath = new OpenLoop();  // forward path boxes
    private OpenLoop closedPath = new OpenLoop();   // full closed path boxes

    private ArrayList<BlackBox> feedbackPath = new ArrayList<BlackBox>(); // feedback path boxes
    private int nFeedbackBoxes = 0;             // number of boxes in feedback path

    private boolean checkPath = false;          // true if segment has been called
    private boolean checkNoMix = true;          // true - no ADC or DAC
    private boolean checkConsolidate = false;   // true if consolidate has been called

    // Constructor
    public ClosedLoop(){
        super("Closed Loop");
    }

    // Add box to the forward path
    public void addBoxToForwardPath(BlackBox box){
        this.forwardPath.addBoxToPath(box);
    }

    // Add box to the open path
    public void addBoxToFeedbackPath(BlackBox box){
        this.feedbackPath.add(box);
        this.nFeedbackBoxes++;
    }

    // Consolidate all boxes into appropriate segments and
    //  combine all boxes into either on forward path box or one closed loop box
    public void consolidate(){

        // add feedback boxes to forward path boxes
        this.closedPath = this.forwardPath.copy();
        for(int i=0; i<this.nFeedbackBoxes; i++){
            this.closedPath.addBoxToPath(this.feedbackPath.get(i));
        }

        // combine forward path boxes
        this.forwardPath.consolidate();

        // combine closed path boxes
        this.closedPath.consolidate();

        // Calculate transfer function
        ComplexPoly fpNumer = this.forwardPath.getSnumer();
        ComplexPoly fpDenom = this.forwardPath.getSdenom();
        ComplexPoly cpNumer = this.closedPath.getSnumer();
        ComplexPoly cpDenom = this.closedPath.getSdenom();
        if(fpDenom.isEqual(cpDenom)){
            super.sNumer = fpNumer.copy();
            this.sDenom = (cpNumer.plus(fpDenom)).copy();
        }
        else{
            super.sNumer = fpNumer.times(cpDenom);
            super.sDenom = cpNumer.plus(cpDenom.times(fpDenom));
        }
        this.checkConsolidate = true;
    }

    // Return number of boxes in the forward path
    public int getNumberOfBoxesInForwardPath(){
        if(!checkConsolidate)this.consolidate();
        return this.forwardPath.getNumberOfBoxes();
    }

    // Return number of boxes in the closed path
    public int getNumberOfBoxesInClosedLoop(){
        if(!checkConsolidate)this.consolidate();
        return this.closedPath.getNumberOfBoxes();
    }

    // Return segment ArrayList for forward path
    public ArrayList<Object> getForwardPathSegmentsArrayList(){
        if(!checkConsolidate)this.consolidate();
        return this.forwardPath.getSegmentsArrayList();
    }

    // Return segment Vector for forward path
    public Vector<Object> getForwardPathSegmentsVector(){
        if(!checkConsolidate)this.consolidate();
        return this.forwardPath.getSegmentsVector();
    }


    // Return segment ArrayList for closed path
    public ArrayList<Object> getClosedLoopSegmentsArrayList(){
        if(!checkConsolidate)this.consolidate();
        return this.closedPath.getSegmentsArrayList();
    }

    // Return segment Vector for closed path
    public Vector<Object> getClosedLoopSegmentsVector(){
        if(!checkConsolidate)this.consolidate();
        return this.closedPath.getSegmentsVector();
    }

   // Return number of segments in the forward path
    public int getNumberOfSegmentsInForwardPath(){
        if(!checkConsolidate)this.consolidate();
        return this.forwardPath.getNumberOfSegments();
    }

    // Return number of segments in the closed path
    public int getNumberOfSegmentsInClosedLoop(){
        if(!checkConsolidate)this.consolidate();
        return this.closedPath.getNumberOfSegments();
    }

    // Return name of all boxes in forward path
    public String getNamesOfBoxesInForwardPath(){
        if(!checkConsolidate)this.consolidate();
        return this.forwardPath.getNamesOfBoxes();
    }

    // Return name of all boxes in closed path
    public String getNamesOfBoxesInClosedLoop(){
        if(!checkConsolidate)this.consolidate();
        return this.closedPath.getNamesOfBoxes();
    }

    // Remove all boxes from the path
    public void removeAllBoxes(){
        this.forwardPath.removeAllBoxes();
        this.closedPath.removeAllBoxes();
    }

    // Deep copy
    public ClosedLoop copy(){
        if(this==null){
            return null;
        }
        else{
            ClosedLoop bb = new ClosedLoop();

            bb.nFeedbackBoxes = this.nFeedbackBoxes;
            bb.checkPath = this.checkPath;
            bb.checkNoMix = this.checkNoMix;
            bb.checkConsolidate = this.checkConsolidate;
            bb.forwardPath = this.forwardPath.copy();
            bb.closedPath = this.closedPath.copy();
            if(this.feedbackPath.size()==0){
                bb.feedbackPath = new ArrayList<BlackBox>();
            }
            else{
                for(int i=0; i<feedbackPath.size(); i++)bb.feedbackPath.add((this.feedbackPath.get(i)).copy());
            }

            bb.sampLen = this.sampLen;
            bb.inputT = this.inputT.clone();
            bb.outputT = this.outputT.clone();
            bb.time = this.time.clone();
            bb.forgetFactor = this.forgetFactor;
            bb.deltaT = this.deltaT;
            bb.sampFreq = this.sampFreq;
            bb.inputS = this.inputS.copy();
            bb.outputS = this.outputS.copy();
            bb.sValue = this.sValue.copy();
            bb.zValue = this.zValue.copy();
            bb.sNumer = this.sNumer.copy();
            bb.sDenom = this.sDenom.copy();
            bb.zNumer = this.zNumer.copy();
            bb.zDenom = this.zDenom.copy();
            bb.sPoles = Complex.copy(this.sPoles);
            bb.sZeros = Complex.copy(this.sZeros);
            bb.zPoles = Complex.copy(this.zPoles);
            bb.zZeros = Complex.copy(this.zZeros);
            bb.sNumerDeg = this.sNumerDeg;
            bb.sDenomDeg = this.sDenomDeg;
            bb.zNumerDeg = this.zNumerDeg;
            bb.zDenomDeg = this.zDenomDeg;
            bb.deadTime = this.deadTime;
            bb.orderPade = this.orderPade;
            bb.sNumerPade = this.sNumerPade.copy();
            bb.sDenomPade = this.sDenomPade.copy();
            bb.sPolesPade = Complex.copy(this.sPolesPade);
            bb.sZerosPade = Complex.copy(this.sZerosPade);
            bb.sNumerDegPade = this.sNumerDegPade;
            bb.sDenomDegPade = this.sDenomDegPade;
            bb.maptozero = this.maptozero;
            bb.padeAdded = this.padeAdded;
            bb.integrationSum = this.integrationSum;
            bb.integMethod = this.integMethod;
            bb.ztransMethod = this.ztransMethod;
            bb.name = this.name;
            bb.fixedName = this.fixedName;
            bb.nPlotPoints = this.nPlotPoints;

            return bb;
        }
    }

    // Clone - overrides Java.Object method clone
    public Object clone(){
        Object ret = null;

        if(this!=null){
            ClosedLoop bb = new ClosedLoop();

            bb.nFeedbackBoxes = this.nFeedbackBoxes;
            bb.checkPath = this.checkPath;
            bb.checkNoMix = this.checkNoMix;
            bb.checkConsolidate = this.checkConsolidate;
            bb.forwardPath = this.forwardPath.copy();
            bb.closedPath = this.closedPath.copy();
            if(this.feedbackPath.size()==0){
                bb.feedbackPath = new ArrayList<BlackBox>();
            }
            else{
                for(int i=0; i<feedbackPath.size(); i++)bb.feedbackPath.add((this.feedbackPath.get(i)).copy());
            }

            bb.sampLen = this.sampLen;
            bb.inputT = this.inputT.clone();
            bb.outputT = this.outputT.clone();
            bb.time = this.time.clone();
            bb.forgetFactor = this.forgetFactor;
            bb.deltaT = this.deltaT;
            bb.sampFreq = this.sampFreq;
            bb.inputS = this.inputS.copy();
            bb.outputS = this.outputS.copy();
            bb.sValue = this.sValue.copy();
            bb.zValue = this.zValue.copy();
            bb.sNumer = this.sNumer.copy();
            bb.sDenom = this.sDenom.copy();
            bb.zNumer = this.zNumer.copy();
            bb.zDenom = this.zDenom.copy();
            bb.sPoles = Complex.copy(this.sPoles);
            bb.sZeros = Complex.copy(this.sZeros);
            bb.zPoles = Complex.copy(this.zPoles);
            bb.zZeros = Complex.copy(this.zZeros);
            bb.sNumerDeg = this.sNumerDeg;
            bb.sDenomDeg = this.sDenomDeg;
            bb.zNumerDeg = this.zNumerDeg;
            bb.zDenomDeg = this.zDenomDeg;
            bb.deadTime = this.deadTime;
            bb.orderPade = this.orderPade;
            bb.sNumerPade = this.sNumerPade.copy();
            bb.sDenomPade = this.sDenomPade.copy();
            bb.sPolesPade = Complex.copy(this.sPolesPade);
            bb.sZerosPade = Complex.copy(this.sZerosPade);
            bb.sNumerDegPade = this.sNumerDegPade;
            bb.sDenomDegPade = this.sDenomDegPade;
            bb.maptozero = this.maptozero;
            bb.padeAdded = this.padeAdded;
            bb.integrationSum = this.integrationSum;
            bb.integMethod = this.integMethod;
            bb.ztransMethod = this.ztransMethod;
            bb.name = this.name;
            bb.fixedName = this.fixedName;
            bb.nPlotPoints = this.nPlotPoints;

            ret = (Object) bb;
        }
        return ret;
    }

}



