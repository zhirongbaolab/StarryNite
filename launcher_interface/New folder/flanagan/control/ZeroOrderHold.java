/*          Class ZeroOrderHold
*
*           This class contains the constructor to create an instance of
*           a zero order hold (ZOH) and the methods needed to use this ZOH
*           in control loops in the time domain, Laplace transform s domain
*           or the z-transform z domain.
*
*           s-domain transfer function = (1 - exp(-Td.s))/s
*           Td is the delay time.
*           Pade approximation always used in s-domain
*           1 to 4 order Pade approximations available
*
*           This class is a subclass of the superclass BlackBox.
*
*           Author:  Michael Thomas Flanagan.
*
*           Created: 26 June 2003.
*           Updated: 2 July 2006, 6 April 2008, 2 December 2008
*
*           DOCUMENTATION:
*           See Michael T Flanagan's JAVA library on-line web page:
*           http://www.ee.ucl.ac.uk/~mflanaga/java/ZeroOrderHold.html
*           http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*           Copyright (c) 2003 - 2008 Michael Thomas Flanagan
*
*           PERMISSION TO COPY:
*
*           Permission to use, copy and modify this software and its documentation for NON-COMMERCIAL purposes is
*           provided that an acknowledgement to the author, Dr Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies
*           and associated documentation or publications.
*
*           Redistributions of the source code of this source code, or parts of the source codes, must retain the above copyright notice,
*           this list of conditions and the following disclaimer and requires written permission from the Michael Thomas Flanagan:
*
*           Redistribution in binary form of all or parts of this class must reproduce the above copyright notice, this list of conditions and
*           the following disclaimer in the documentation and/or other materials provided with the distribution and requires written permission
*           from the Michael Thomas Flanagan:
*
*           Dr Michael Thomas Flanagan makes no representations about the suitability or fitness of the software for any or for a particular purpose.
*           Dr Michael Thomas Flanagan shall not be liable for any damages suffered as a result of using, modifying or distributing this software
*           or its derivatives.
*
***************************************************************************************/

package flanagan.control;

import flanagan.complex.Complex;
import flanagan.complex.ComplexPoly;

public class ZeroOrderHold extends BlackBox{

    // Constructor
    public ZeroOrderHold(double deltaT, int orderPade){
        super("ZeroOrderHold");
        super.sPoles = Complex.oneDarray(1);
        super.setDeltaT(deltaT);
        super.setPadeOrder(orderPade);
        this.setNumDen(deltaT);

    }

    // Constructor
    // Default Pade approximation order = 2
    public ZeroOrderHold(double deltaT){
        super("ZeroOrderHold");
        super.sPoles = Complex.oneDarray(1);
        super.setDeltaT(deltaT);
        this.setNumDen(deltaT);
    }

    // set the numerators and denominators
    // same polynomials, using Pade approximation for Pade and non-Pade forms
    public void setNumDen(double deltaT){
        // set denominator, s
        super.sDenom = new ComplexPoly(0.0D, 1.0D);
        super.sPoles[0].reset(0.0D, 0.0D);

        // set exp(-sT) part of pade numerator
        super.sNumer = new ComplexPoly(1.0D);
        super.deadTime = deltaT;
        super.pade();
        super.deadTime = 0.0D;

        // add 1 to exp(-sT)[=padeNumer/padeDenom]/s
        super.sNumerPade = super.sNumerPade.plus(super.sDenomPade);
        super.sZerosPade = sNumerPade.roots();
    }

    // Deep copy
    public ZeroOrderHold copy(){
        if(this==null){
            return null;
        }
        else{
            ZeroOrderHold bb = new ZeroOrderHold(this.deltaT, this.orderPade);

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
            ZeroOrderHold bb = new ZeroOrderHold(this.deltaT, this.orderPade);

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






