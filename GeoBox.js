
//Declare the Constants for Calculations
var Coords  = function(Datum) {
	if(!Datum){
		alert('No Datum Selected. You Must Specify a Datum First');
		return false;
	}
	this.Datum = Datum;
};


Coords.prototype.MakeDigraph = function(y, utmz) {
	//Inputs y utmz
    //alert(utmz);
    Letr = Math.floor((this.utmz-1)*8 + (x)/100000);
    Letr = Letr - 24*Math.floor(Letr/24)-1;
    Digraph = this.DigraphLetrsE.charAt(Letr);
    //alert("x=   "+x);
    //alert(DigraphLetrsE.charAt(Letr));
    //First (Easting) Character Found
    Letr = Math.floor(y/100000);
    //Odd zones start with A at equator, even zones with F
    if (utmz/2 == Math.floor(utmz/2)){Letr = Letr+5;}
    Letr = Letr - 20*Math.floor(Letr/20);
    Digraph = Digraph + DigraphLetrsN.charAt(Letr);
};

Coords.prototype.DDtoDMS = function(yd,xd) {
    //Input= xd(long) and yd(lat)
    //Output = xdd xm xs (long) and ydd ym ys (lat)
    var ydd = Math.floor(Math.abs(yd));
    var ym = Math.floor(60*(Math.abs(yd) - ydd));
    var ys = 3600*(Math.abs(yd)-ydd - ym/60);
    if (yd<0){ydd=-ydd;}
    var xdd = Math.floor(Math.abs(xd));
    var xm = Math.floor(60*(Math.abs(xd) - xdd));
    var xs = 3600*(Math.abs(xd)-xdd - xm/60);
    if (xd<0){xdd=-xdd;}
    var DMS = {
	  	LatDD: ydd,
	  	LatDM: ym,
	  	LatDS: ys,
	  	LongDD: xdd,
	  	LongDM: xm,
	  	LongDS: xs
	};
	return DMS;
};//End DDtoDMS

Coords.prototype.DMStoDD = function(ydd,ym,ys,xdd,xm,xs) {
	//Input = xdd xm xs (long) and ydd ym ys (lat)
    //Output= xd(long) and yd(lat)
    var xd = Math.abs(xdd) + xm/60 + xs/3600;
    var yd = Math.abs(ydd) + ym/60 + ys/3600;
    if (ydd<0){yd=-yd;}
    if (xdd<0){xd=-xd;}
    var DD = {
	  	Lat: yd,
	  	Long: xd,
	};
	return DD;
};//End DMStoDD

Coords.prototype.DDtoDegM = function(yd,xd) {
    //Input= xd(long) and yd(lat)
    //Output = xdd xdm  (long) and ydd ydm (lat)
    var ydd = Math.floor(Math.abs(yd));
    var ydm = (60*(Math.abs(yd) - ydd));
    if (yd<0){ydd=-ydd;}
    var xdd = Math.floor(Math.abs(xd));
    var xdm = (60*(Math.abs(xd) - xdd));
    if (xd<0){xdd=-xdd;}
    var DegM = {
	  	LatDD: ydd,
	  	LatDM: ydm,
	  	LongDD: xdd,
	  	LongDM: xdm,
	};
	return DegM;
};//End DDtoDegM

Coords.prototype.DMStoDegM = function(ydd,ym,ys,xdd,xm,xs) {
	//Input = xdd xm xs (long) and ydd ym ys (lat)
    //Output= xd(long) and yd(lat)
    var	xdm = xs/60 + xm;
    var	ydm = ys/60 + ym;
 	var DegM = {
	  	LatDD: ydd,
	  	LatDM: ydm,
	  	LongDD: xdd,
	  	LongDM: xdm,
	};
	return DegM;
};//End DMStoDegM

Coords.prototype.DegMtoDD = function(ydd,ydm,xdd,xdm) {
	//Input = xdd xdm (long) and ydd ydm (lat)
    //Output= xd(long) and yd(lat)
    var xd = xdd + xdm/60;
    var yd = ydd + ydm/60;
    var DD = {
	  	Lat: yd,
	  	Long: xd,
	};
	return DD;
};//End DegMtoDD

Coords.prototype.DegMtoDMS = function(ydd,ydm,xdd,xdm) {
	//Input = xdd xdm (long) and ydd ydm (lat)
    //Output= xd(long) and yd(lat)
    var xm = Math.floor(Math.abs(xdm));
    var xs = (60*(Math.abs(xdm) - xm));
    var ym = Math.floor(Math.abs(ydm));
    var ys = (60*(Math.abs(ydm) - ym));
	var DMS = {
	  	LatDD: ydd,
	  	LatDM: ym,
	  	LatDS: ys,
	  	LongDD: xdd,
	  	LongDM: xm,
	  	LongDS: xs
	};
	return DMS;
};//End DegMtoDMS

Coords.prototype.DDtoUTM = function(latd,lngd) {
	return this._GeoToUTM(latd,lngd);
};

Coords.prototype.DMStoUTM = function(ydd,ym,ys,xdd,xm,xs) {
	var DD = this.DMStoDD(ydd,ym,ys,xdd,xm,xs);
	return this._GeoToUTM(DD['Lat'],DD['Long']);
};
Coords.prototype.DegMtoUTM = function(ydd,ydm,xdd,xdm) {
	var DD = this.DegMtoDD(ydd,ydm,xdd,xdm);
	return this._GeoToUTM(DD['Lat'],DD['Long']);
};

Coords.prototype._GeoToUTM = function(latd,lngd) {
	var latd = parseFloat(latd);
	var lngd = parseFloat(lngd);
	var Datum = this.Datum;
    //Convert Latitude and Longitude to UTM
    /*
    Declarations();
	Available Datums
		WGS_84
    	NAD_83
    	GRS_80
    	WGS_72
    	Australian_1965
    	Krasovsky_1940
    	North_American_1927
    	International_1924
    	Hayford_1909
    	Clarke_1880
    	Clarke_1866
    	Airy_1830
    	Bessel_1841
    	Everest_1830
	*/

    //Symbols as used in USGS PP 1395: Map Projections - A Working Manual
    var DatumEqRad = {
    	WGS_84: 6378137.0,
    	NAD_83: 6378137.0,
    	GRS_80: 6378137.0,
    	WGS_72: 6378135.0,
    	Australian_1965: 6378160.0,
    	Krasovsky_1940: 6378245.0,
    	North_American_1927: 6378206.4,
    	International_1924: 6378388.0,
    	Hayford_1909: 6378388.0,
    	Clarke_1880: 6378249.1,
    	Clarke_1866: 6378206.4,
    	Airy_1830: 6377563.4,
    	Bessel_1841: 6377397.2,
    	Everest_1830: 6377276.3
    }
    var DatumFlat = {
    	WGS_84: 298.2572236, 
		NAD_83: 298.2572236, 
		GRS_80: 298.2572215, 
		WGS_72: 298.2597208, 
		Australian_1965: 298.2497323, 
		Krasovsky_1940: 298.2997381, 
		North_American_1927: 294.9786982,
	    International_1924: 296.9993621,
	    Hayford_1909: 296.9993621, 
	    Clarke_1880: 293.4660167, 
	    Clarke_1866: 294.9786982, 
	    Airy_1830: 299.3247788, 
	    Bessel_1841: 299.1527052, 
	    Everest_1830: 300.8021499
    }

    var k0 = 0.9996;//scale on central meridian
    var a = parseFloat(DatumEqRad[Datum]);//equatorial radius, meters. 
    var f = 1/parseFloat(DatumFlat[Datum]);//polar flattening.
    var b = a*(1-f);//polar axis.
  
    var e = Math.sqrt(1 - b*b/a*a);//eccentricity
    var drad = Math.PI/180;//Convert degrees to radians)
    var phi = 0;//latitude (north +, south -), but uses phi in reference
    var e0 = e/Math.sqrt(1 - e*e);//e prime in reference
    var N = a/Math.sqrt(1-Math.pow(e*Math.sin(phi)),2);
    var T = Math.pow(Math.tan(phi),2);
    var C = Math.pow(e*Math.cos(phi),2);
    var lng = 0;//Longitude (e = +, w = -) - can't use long - reserved word
    var lng0 = 0;//longitude of central meridian
    var M = 0;//M requires calculation
    var x = 0;//x coordinate
    var y = 0;//y coordinate
    var k = 1;//local scale
    var utmz = 30;//utm zone
    var zcm = 0;//zone central meridian
    var DigraphLetrsE = "ABCDEFGHJKLMNPQRSTUVWXYZ";
    var DigraphLetrsN = "ABCDEFGHJKLMNPQRSTUV";
    var OOZok = false;

    var e = Math.sqrt(1 - (b/a)*(b/a));//eccentricity    
  	console.log(a);
    console.log(f);
    console.log(b);
    console.log(e);
    if(isNaN(latd)|| isNaN(lngd)){
        alert("Non-Numeric Input Value");
    }
    if(latd <-90 || latd> 90){
        alert("Latitude must be between -90 and 90");
    }
    if(lngd <-180 || lngd > 180){
        alert("Latitude must be between -180 and 180");
    }
 
    var xd = lngd;
    var yd = latd;
    
    var phi = latd*drad;//Convert latitude to radians
    var lng = lngd*drad;//Convert longitude to radians
    var utmz = 1 + Math.floor((lngd+180)/6);//calculate utm zone
    var latz = 0;//Latitude zone: A-B S of -80, C-W -80 to +72, X 72-84, Y,Z N of 84
    if (latd > -80 && latd < 72){latz = Math.floor((latd + 80)/8)+2;}
    if (latd > 72 && latd < 84){latz = 21;}
    if (latd > 84){latz = 23;}
        
    var zcm = 3 + 6*(utmz-1) - 180;//Central meridian of zone
    //Calculate Intermediate Terms
    var e0 = e/Math.sqrt(1 - e*e);//Called e prime in reference
    var esq = (1 - (b/a)*(b/a));//e squared for use in expansions
    var e0sq = e*e/(1-e*e);// e0 squared - always even powers
    var N = a/Math.sqrt(1-Math.pow(e*Math.sin(phi),2));
    var T = Math.pow(Math.tan(phi),2);
    var C = e0sq*Math.pow(Math.cos(phi),2);
    var A = (lngd-zcm)*drad*Math.cos(phi);
    //Calculate M
    var M = phi*(1 - esq*(1/4 + esq*(3/64 + 5*esq/256)));
    M = M - Math.sin(2*phi)*(esq*(3/8 + esq*(3/32 + 45*esq/1024)));
    M = M + Math.sin(4*phi)*(esq*esq*(15/256 + esq*45/1024));
    M = M - Math.sin(6*phi)*(esq*esq*esq*(35/3072));
    M = M*a;//Arc length along standard meridian
    var M0 = 0;//M0 is M for some origin latitude other than zero. Not needed for standard UTM

    //Calculate UTM Values
    var x = k0*N*A*(1 + A*A*((1-T+C)/6 + A*A*(5 - 18*T + T*T + 72*C -58*e0sq)/120));//Easting relative to CM
   	x=x+500000;//Easting standard 
    var y = k0*(M - M0 + N*Math.tan(phi)*(A*A*(1/2 + A*A*((5 - T + 9*C + 4*C*C)/24 + A*A*(61 - 58*T + T*T + 600*C - 330*e0sq)/720))));//Northing from equator
    var yg = y + 10000000;//yg = y global, from S. Pole
    if (y < 0){y = 10000000+y;}
    var SHem = (phi<0) ? true : false;
    //Output into UTM Boxes
    var UTM = {
	  	UTMz: utmz,
	  	UTMe: Math.round(10*(x))/10,
	  	UTMn: Math.round(10*y)/10,
	  	SHem: SHem,
	};
    /*//NATO UTM
    document.getElementById("UTMLonZoneBox2").value = utmz;
    document.getElementById("UTMLatZoneBox2").value = DigraphLetrsE[latz];
    document.getElementById("UTMeBox2").value = Math.round(10*(x-100000*Math.floor(x/100000)))/10;
    document.getElementById("UTMnBox2").value = Math.round(10*(y-100000*Math.floor(y/100000)))/10;
//Generate Digraph
    MakeDigraph();
    document.getElementById("UTMDgBox2").value = Digraph;*/
    return UTM;
}//close Geog to UTM
///////////////////////////////////////////////////////////////////////


Coords.prototype.UTMtoDD = function(UTMe,UTMn,UTMz,SHem) {
    return this._UTMtoGeog(UTMe,UTMn,UTMz,SHem);
};
Coords.prototype.UTMtoDMS = function(UTMe,UTMn,UTMz,SHem) {
    var DD = this._UTMtoGeog(UTMe,UTMn,UTMz,SHem);
    return this.DDtoDMS(DD['Lat'],DD['Long']); 
};
Coords.prototype.UTMtoDegM = function(UTMe,UTMn,UTMz,SHem) {
    var DD = this._UTMtoGeog(UTMe,UTMn,UTMz,SHem);
    return this.DDtoDegM(DD['Lat'],DD['Long']); 
};

Coords.prototype._UTMtoGeog = function(UTMe,UTMn,UTMz,SHem) {
    //Convert UTM Coordinates to Geographic
    var Datum = this.Datum;
    //Convert Latitude and Longitude to UTM
    /*
    Declarations();
    Available Datums
        WGS_84
        NAD_83
        GRS_80
        WGS_72
        Australian_1965
        Krasovsky_1940
        North_American_1927
        International_1924
        Hayford_1909
        Clarke_1880
        Clarke_1866
        Airy_1830
        Bessel_1841
        Everest_1830
    */

    //Symbols as used in USGS PP 1395: Map Projections - A Working Manual
    var DatumEqRad = {
        WGS_84: 6378137.0,
        NAD_83: 6378137.0,
        GRS_80: 6378137.0,
        WGS_72: 6378135.0,
        Australian_1965: 6378160.0,
        Krasovsky_1940: 6378245.0,
        North_American_1927: 6378206.4,
        International_1924: 6378388.0,
        Hayford_1909: 6378388.0,
        Clarke_1880: 6378249.1,
        Clarke_1866: 6378206.4,
        Airy_1830: 6377563.4,
        Bessel_1841: 6377397.2,
        Everest_1830: 6377276.3
    }
    var DatumFlat = {
        WGS_84: 298.2572236, 
        NAD_83: 298.2572236, 
        GRS_80: 298.2572215, 
        WGS_72: 298.2597208, 
        Australian_1965: 298.2497323, 
        Krasovsky_1940: 298.2997381, 
        North_American_1927: 294.9786982,
        International_1924: 296.9993621,
        Hayford_1909: 296.9993621, 
        Clarke_1880: 293.4660167, 
        Clarke_1866: 294.9786982, 
        Airy_1830: 299.3247788, 
        Bessel_1841: 299.1527052, 
        Everest_1830: 300.8021499
    }

    var k0 = 0.9996;//scale on central meridian
    var a = parseFloat(DatumEqRad[Datum]);//equatorial radius, meters. 
    var f = 1/parseFloat(DatumFlat[Datum]);//polar flattening.
    var b = a*(1-f);//polar axis.
  
    var e = Math.sqrt(1 - b*b/a*a);//eccentricity
    var drad = Math.PI/180;//Convert degrees to radians)
    var phi = 0;//latitude (north +, south -), but uses phi in reference
    var e0 = e/Math.sqrt(1 - e*e);//e prime in reference
    var N = a/Math.sqrt(1-Math.pow(e*Math.sin(phi)),2);
    var T = Math.pow(Math.tan(phi),2);
    var C = Math.pow(e*Math.cos(phi),2);
    var lng = 0;//Longitude (e = +, w = -) - can't use long - reserved word
    var lng0 = 0;//longitude of central meridian
    var M = 0;//M requires calculation
    var x = 0;//x coordinate
    var y = 0;//y coordinate
    var k = 1;//local scale
    var utmz = 30;//utm zone
    var zcm = 0;//zone central meridian
    var DigraphLetrsE = "ABCDEFGHJKLMNPQRSTUVWXYZ";
    var DigraphLetrsN = "ABCDEFGHJKLMNPQRSTUV";
    var OOZok = false;


    
    k0 = 0.9996;//scale on central meridian
    b = a*(1-f);//polar axis.
    e = Math.sqrt(1 - (b/a)*(b/a));//eccentricity
    e0 = e/Math.sqrt(1 - e*e);//Called e prime in reference
    esq = (1 - (b/a)*(b/a));//e squared for use in expansions
    e0sq = e*e/(1-e*e);// e0 squared - always even powers
    x = parseFloat(UTMe);
    if (x<160000 || x>840000){alert("Outside permissible range of easting values \n Results may be unreliable \n Use with caution");} 
    y = parseFloat(UTMn);
    if (y<0){alert("Negative values not allowed \n Results may be unreliable \n Use with caution");}
    if (y>10000000){alert("Northing may not exceed 10,000,000 \n Results may be unreliable \n Use with caution");}
    utmz = parseFloat(UTMz);
    zcm = 3 + 6*(utmz-1) - 180;//Central meridian of zone
    var e1 = (1 - Math.sqrt(1 - e*e))/(1 + Math.sqrt(1 - e*e));//Called e1 in USGS PP 1395 also
    var M0 = 0;//In case origin other than zero lat - not needed for standard UTM
    M = M0 + y/k0;//Arc length along standard meridian. 
    if (SHem == true){M=M0+(y-10000000)/k;}
    var mu = M/(a*(1 - esq*(1/4 + esq*(3/64 + 5*esq/256))));
    var phi1 = mu + e1*(3/2 - 27*e1*e1/32)*Math.sin(2*mu) + e1*e1*(21/16 -55*e1*e1/32)*Math.sin(4*mu);//Footprint Latitude
    phi1 = phi1 + e1*e1*e1*(Math.sin(6*mu)*151/96 + e1*Math.sin(8*mu)*1097/512);

    var C1 = e0sq*Math.pow(Math.cos(phi1),2);
    var T1 = Math.pow(Math.tan(phi1),2);
    var N1 = a/Math.sqrt(1-Math.pow(e*Math.sin(phi1),2));
    var R1 = N1*(1-e*e)/(1-Math.pow(e*Math.sin(phi1),2));
    var D = (x-500000)/(N1*k0);
    var phi = (D*D)*(1/2 - D*D*(5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*e0sq)/24);
    phi = phi + Math.pow(D,6)*(61 + 90*T1 + 298*C1 + 45*T1*T1 -252*e0sq - 3*C1*C1)/720;
    phi = phi1 - (N1*Math.tan(phi1)/R1)*phi;
    
//Output Latitude
    var latitude = Math.floor(1000000*phi/drad)/1000000;
        
//Longitude
    var lng = D*(1 + D*D*((-1 -2*T1 -C1)/6 + D*D*(5 - 2*C1 + 28*T1 - 3*C1*C1 +8*e0sq + 24*T1*T1)/120))/Math.cos(phi1);
    var lngd = zcm+lng/drad;
    
//Output Longitude
    var longitude = Math.floor(1000000*lngd)/1000000;

    var DD = {
        Lat: latitude,
        Long: longitude,
    };

    return DD;

    };