/*      Class OpenLoop
*
*       This class supports the creation of a path of Black Boxes
*       i.e. of instances of BlackBox and of any of its subclasses,
*       e.g. PropIntDeriv, FirstOrder, and the methods to combine
*       these into both a single instance of BlackBox and a Vector
*       of analogue segments, digital segments and converters.
*
*       Author:  Michael Thomas Flanagan.
*
*       Created: August 2002
*	    Updated: 12 July 2003, 10 May 2005, 2 July 2006, 7 June 2007, 6 April 2008, 5 July 2008
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

import java.util.ArrayList;
import java.util.Vector;
import flanagan.complex.Complex;
import flanagan.complex.ComplexPoly;

public class OpenLoop extends BlackBox{
    private ArrayList<BlackBox> openPath = new ArrayList<BlackBox>(); // open path boxes
    private ArrayList<Object> segments = new ArrayList<Object>();     // start of segment, end of segment, type of each segment, i.e. analogue, digital, AtoD, DtoA, ZOH

    private int nBoxes = 0;                     // number of boxes in original path
    private int nSeg = 0;                       // number of analogue, digital, AtoD, ZOH segments

    private boolean checkPath = false;          // true if segment has been called
    private boolean checkNoMix = true;          // true - no ADC or DAC
    private boolean checkConsolidate = false;   // true if consolidate has been called

    // Constructor
    public OpenLoop(){
        super("OpenLoop");
    }

    // Add box to the open path
    public void addBoxToPath(BlackBox box){
        this.openPath.add(box);
        this.nBoxes++;
    }

    // Consolidate all boxes into appropriate segments and combine all boxes into one box
    public void consolidate(){
        // Empty segments ArrayList if openPath ArrayList has been updated
        if(!segments.isEmpty()){
            segments.clear();
            this.nBoxes = 0;
            this.nSeg = 0;
            this.checkNoMix = true;
            this.checkPath = false;
        }

        // Find analogue, digital and conversion segments in OpenLoop
        this.segment();
        // Combine all boxes into a single box and make this instance that combined box
        BlackBox aa = null;
        if(this.nSeg==1){
            aa = (BlackBox) this.segments.get(3);
        }
        else{
            aa = this.combineSegment(0, this.nBoxes);
        }
        super.sNumer = aa.sNumer.copy();
        super.sDenom = aa.sDenom.copy();
        super.sNumerPade = aa.sNumerPade.copy();
        super.sDenomPade = aa.sDenomPade.copy();
        super.sNumerDeg = aa.sNumerDeg;
        super.sDenomDeg = aa.sDenomDeg;
        super.sNumerDegPade = aa.sNumerDegPade;
        super.sDenomDegPade = aa.sDenomDegPade;
        super.deadTime = aa.deadTime;
        super.sZeros = Complex.copy(aa.sZeros);
        super.sPoles = Complex.copy(aa.sPoles);
        super.sZerosPade = Complex.copy(aa.sZerosPade);
        super.sPolesPade = Complex.copy(aa.sPolesPade);
        super.padeAdded=true;

        this.checkConsolidate = true;
    }

    // Find analogue and digital segments
    public void segment(){
        this.checkPath = true;              // this method, segment, has been called
        this.nBoxes = openPath.size();      // number of boxes in openPath
        String thisName = " ";              // name of current openPath box under examination

        // Find analogue, digital, ZOH/ADC and DAC segments
        int iStart1=0;                      // start index of block under examination
        int iEnd1=0;                        // final index of block under examination
        int iStart2=0;                      // start index of next block after block under examination
                                            //  to accomodate separate block for ZOH/ADC or DAC conversion
        int iEnd2=0;                        // final index of next block after block under examination
                                            //  to accomodate separate block for ZOH/ADC or DAC conversion
        int iNewStart=0;                    // start index of next block after segment/s completed
        int nnBoxes = 0;                    // number of segments in a completed segement test
        int nInSeg = 0;                     // counter for number of openPath boxes in a segment

        String name1 = " ";                 // name of first segment resulting from a segment test
        String name2 = " ";                 // name of second segment resulting from a segment test
        String lastConv = " ";              // name of the last convertor, e.g. AtoD or DtoA, found

        int ii = 0;                         // counter indicating position along openPath

        double deltaThold = 0.0D;           // holds value of deltaT for box under examination

        while(ii<this.nBoxes){
            nInSeg++;
            BlackBox bb = openPath.get(ii);
            thisName = bb.fixedName;
            //if(bb.deltaT!=0.0D)deltaThold=bb.deltaT;

            // Look for ZOH
            if(thisName.equals("ZeroOrderHold")){
                if(!lastConv.equals(" "))this.checkNoMix=false;
                if(ii<this.nBoxes-1){
                    BlackBox cc = openPath.get(ii+1);

                    // Look for following ADC
                    if(cc.fixedName.equals("AtoD")){
                        if(lastConv.equals("AtoD"))throw new IllegalArgumentException("Two consecutive ADCs with no intervening DAC");
                        if(nInSeg>1){
                            iEnd1 = ii-1;
                            name1 = "analogue";
                            iStart2 = ii;
                            iEnd2 = ii+1;
                            name2 = "AtoD";
                            nnBoxes = 2;
                            this.nSeg += 2;
                            ii = ii+2;
                            iNewStart = iEnd2+1;
                        }
                        else{
                            iEnd1 = ii+1;
                            name1 = "AtoD";
                            nnBoxes = 1;
                            this.nSeg += 1;
                            ii = ii+2;
                            iNewStart = iEnd1+1;
                        }
                        lastConv = "AtoD";
                        nInSeg=0;
                    }
                    else{
                        System.out.println("WARNING!! OpenLoop.checkPath: ZOH without a following ADC");
                        if(nInSeg>1){
                            iEnd1=ii-1;
                            name1 = "analogue";
                            iStart2 = ii;
                            iEnd2 = ii;
                            name2 = "ZOH";
                            nnBoxes = 2;
                            this.nSeg =+ 2;
                            ii = ii+1;
                            iNewStart = iEnd2+1;
                        }
                        else{
                            iEnd1=ii;
                            name1 = "ZOH";
                            nnBoxes = 1;
                            this.nSeg =+ 1;
                            ii = ii+1;
                            iNewStart = iEnd1+1;
                        }
                        nInSeg = 0;
                        lastConv = "ZOH";
                    }
                }
                else{
                    System.out.println("WARNING!! OpenLoop.checkPath: path ends with ZOH");
                    if(nInSeg>1){
                        iEnd1=ii-1;
                        name1 = "analogue";
                        iStart2 = ii;
                        iEnd2 = ii;
                        name2 = "ZOH";
                        nnBoxes = 2;
                        this.nSeg += 2;
                        ii = ii+2;
                        iNewStart = iEnd2+1;
                    }
                    else{
                        iEnd1=ii;
                        name1 = "ZOH";
                        nnBoxes = 1;
                        this.nSeg =+ 1;
                        ii = ii+1;
                        iNewStart = iEnd1+1;
                    }
                    lastConv = "ZOH";
                    nInSeg = 0;
                }
            }
            else{
                if(thisName.equals("AtoD")) throw new IllegalArgumentException("ADC without preceeding ZOH");

                // Look for DAC
                if(thisName.equals("DtoA")){
                    if(lastConv.equals("DtoA"))throw new IllegalArgumentException("Two consecutive DACs with no intervening ADC");
                    if(lastConv.equals("ZOH"))throw new IllegalArgumentException("ZOH followed by DAC");
                    if(!lastConv.equals(" "))this.checkNoMix=false;
                    if(nInSeg>1){
                        iEnd1=ii-1;
                        name1 = "digital";
                        iStart2 = 2;
                        iEnd2 = ii;
                        ii = ii+1;
                        iNewStart = iEnd1+1;
                        nnBoxes = 2;
                        this.nSeg =+ 2;
                    }
                    else{
                        iEnd1=ii;
                        name1 = "DtoA";
                        ii = ii+1;
                        iNewStart = iEnd1+1;
                        nnBoxes = 1;
                        this.nSeg =+ 1;
                    }
                    lastConv = "DtoA";
                    nInSeg = 0;
                }
            }

            // Add segment/s found to segments ArrayList as 4 elements:
            // 1. start index, 2. final index, 3. name, 4. all boxes in segment combined as a single box
            if(nnBoxes>0){
                this.segments.add(new Integer(iStart1));
                this.segments.add(new Integer(iEnd1));
                this.segments.add(name1);
                BlackBox dd = this.combineSegment(iStart1, iEnd1);
                this.segments.add(dd);
                if(nnBoxes==2){
                    this.segments.add(new Integer(iStart2));
                    this.segments.add(new Integer(iEnd2));
                    this.segments.add(name2);
                    BlackBox ee = this.combineSegment(iStart2, iEnd2);
                    this.segments.add(ee);
                }
                iStart1=iNewStart;
            }
            else{
                ii++;
            }
            if(ii>=this.nBoxes && ii!=iNewStart){
                iEnd1=ii-1;
                name1 = "analogue";
                if(lastConv.equals("AtoD"))name1 = "digital";
                this.nSeg =+ 1;
                this.segments.add(new Integer(iStart1));
                this.segments.add(new Integer(iEnd1));
                this.segments.add(name1);
                BlackBox ff = this.combineSegment(iStart1, iEnd1);
                this.segments.add(ff);
            }

        }
    }

    // Combine all boxes between iLow and iHigh into one box
    public BlackBox combineSegment(int iLow, int iHigh){
        BlackBox aa = new BlackBox();           // Black Box to be returned
        int nBoxSeg = iHigh - iLow + 1;         // number of boxes in segment
        int[] numDeg = new int[nBoxSeg];        // array of numerator degrees
        int[] denDeg = new int[nBoxSeg];        // array of denominator degrees
        BlackBox bb = openPath.get(iLow);
        if(!bb.padeAdded)bb.transferPolesZeros();
        aa.sNumerPade = bb.sNumerPade.copy();
        aa.sDenomPade = bb.sDenomPade.copy();
        aa.deadTime = bb.deadTime;
        numDeg[0] = bb.sNumerDegPade;
        denDeg[0] = bb.sDenomDegPade;
        aa.sNumerDegPade = numDeg[0];
        aa.sDenomDegPade = denDeg[0];
        for(int i=1; i<nBoxSeg; i++){
            bb = openPath.get(i+iLow);
            if(!bb.padeAdded)bb.transferPolesZeros();
            if(aa.sNumerPade==null){
                if(bb.sNumerPade!=null){
                    aa.sNumerPade = bb.sNumerPade.copy();
                }
            }
            else{
                if(bb.sNumerPade!=null){
                    aa.sNumerPade = bb.sNumerPade.times(aa.sNumerPade);
                }
            }

            if(aa.sDenomPade==null){
                if(bb.sDenomPade!=null){
                    aa.sDenomPade = bb.sDenomPade.copy();
                }
            }
            else{
                if(bb.sDenomPade!=null){
                    aa.sDenomPade = bb.sDenomPade.times(aa.sDenomPade);
                }
            }

            aa.deadTime += bb.deadTime;
            numDeg[i] = bb.sNumerDegPade;
            denDeg[i] = bb.sDenomDegPade;
            aa.sNumerDegPade += numDeg[i];
            aa.sDenomDegPade += denDeg[i];
        }
        if(aa.sNumerDegPade>0){
            aa.sZerosPade = Complex.oneDarray(aa.sNumerDegPade);
            int numK = 0;
            int denK = 0;
            for(int i=0; i<nBoxSeg; i++){
                bb = openPath.get(i+iLow);
                if(bb.sNumerDegPade>0){
                    for(int j=0; j<numDeg[i]; j++){
                        aa.sZerosPade[numK] = bb.sZerosPade[j].copy();
                        numK++;
                    }
                }
            }
        }

        if(aa.sNumerDegPade>0){
            aa.sPolesPade = Complex.oneDarray(aa.sDenomDegPade);
            int numK = 0;
            int denK = 0;
            for(int i=0; i<nBoxSeg; i++){
                bb = openPath.get(i+iLow);
                if(bb.sNumerDegPade>0){
                    for(int j=0; j<denDeg[i]; j++){
                        aa.sPolesPade[denK] = bb.sPolesPade[j].copy();
                        denK++;
                    }
                }
            }
        }

        aa.zeroPoleCancellation();
        aa.padeAdded = true;
        aa.sNumerDeg = aa.sNumerDegPade;
        aa.sDenomDeg = aa.sDenomDegPade;
        aa.sNumer = aa.sNumerPade.copy();
        aa.sNumer = aa.sNumerPade.copy();
        aa.sZeros = Complex.copy(aa.sZerosPade);
        aa.sPoles = Complex.copy(aa.sPolesPade);
        return aa;

    }

    // Return number of boxes in path
    public int getNumberOfBoxes(){
        if(!checkConsolidate)this.consolidate();
        return this.nBoxes;
    }

    // Return segment ArrayList
    public ArrayList<Object> getSegmentsArrayList(){
        if(!checkConsolidate)this.consolidate();
        return this.segments;
    }

    // Return segment Vector
    public Vector<Object> getSegmentsVector(){
        if(!checkConsolidate)this.consolidate();
        ArrayList<Object> seg = this.segments;
        Vector<Object> ret = null;
        if(seg!=null){
            int n = seg.size();
            ret = new Vector<Object>(n);
            for(int i=0; i<n; i++)ret.addElement(seg.get(i));
        }
        return ret;
    }

    // Return number of segments in path
    public int getNumberOfSegments(){
        if(!checkConsolidate)this.consolidate();
        return this.nSeg;
    }

    // Return name of all boxes in path
    public String getNamesOfBoxes(){
        if(!checkConsolidate)this.consolidate();
        String names = "";
        for(int i=0; i<this.nBoxes; i++){
            BlackBox bb = openPath.get(i);
            names = names + i +": "+bb.getName() + "   ";
        }
        return names;
    }

    // Remove all boxes from the path
    public void removeAllBoxes(){
         // Empty openPath ArrayList
        if(!openPath.isEmpty()){
            openPath.clear();
        }

        // Empty segments ArrayList
        if(!segments.isEmpty()){
            segments.clear();
        }
        this.nSeg = 0;
        this.checkNoMix = true;
        this.checkPath = false;
    }

    // Deep copy
    public OpenLoop copy(){
        if(this==null){
            return null;
        }
        else{
            OpenLoop bb = new OpenLoop();

            bb.nBoxes = this.nBoxes;
            bb.nSeg = this.nSeg;
            bb.checkPath = this.checkPath;
            bb.checkNoMix = this.checkNoMix;
            bb.checkConsolidate = this.checkConsolidate;
            if(this.openPath.size()==0){
                bb.openPath = new ArrayList<BlackBox>();
            }
            else{
                for(int i=0; i<openPath.size(); i++)bb.openPath.add((this.openPath.get(i)).copy());
            }
            if(this.segments.size()==0){
                bb.segments = new ArrayList<Object>();
            }
            else{
                int j=0;
                for(int i=0; i<this.nSeg; i++){
                    Integer holdI1 = (Integer)this.segments.get(j);
                    int ii = holdI1.intValue();
                    bb.segments.add(new Integer(ii));
                    j++;
                    Integer holdI2 = (Integer)this.segments.get(j);
                    ii = holdI2.intValue();
                    bb.segments.add(new Integer(ii));
                    j++;
                    String holdS = (String)this.segments.get(j);
                    bb.segments.add(holdS);
                    j++;
                    bb.segments.add(((BlackBox)this.segments.get(j)).copy());
                    j++;
                }
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
            OpenLoop bb = new OpenLoop();

            bb.nBoxes = this.nBoxes;
            bb.nSeg = this.nSeg;
            bb.checkPath = this.checkPath;
            bb.checkNoMix = this.checkNoMix;
            bb.checkConsolidate = this.checkConsolidate;
            if(this.openPath.size()==0){
                bb.openPath = new ArrayList<BlackBox>();
            }
            else{
                for(int i=0; i<openPath.size(); i++)bb.openPath.add((this.openPath.get(i)).copy());
            }
            if(this.segments.size()==0){
                bb.segments = new ArrayList<Object>();
            }
            else{
                int j=0;
                for(int i=0; i<this.nSeg; i++){
                    Integer holdI1 = (Integer)this.segments.get(j);
                    int ii = holdI1.intValue();
                    bb.segments.add(new Integer(ii));
                    j++;
                    Integer holdI2 = (Integer)this.segments.get(j);
                    ii = holdI2.intValue();
                    bb.segments.add(new Integer(ii));
                    j++;
                    String holdS = (String)this.segments.get(j);
                    bb.segments.add(holdS);
                    j++;
                    bb.segments.add(((BlackBox)this.segments.get(j)).copy());
                    j++;
                }
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


