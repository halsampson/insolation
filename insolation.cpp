// insolation.cpp

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <windows.h>

const double PI = 3.14159265359;
#define radians(degrees)  (PI * (degrees) / 180)

#include "private.txt"  // #define LatitudeDegrees and LongitudeDegrees

const double latitude = radians(LatitudeDegrees);
const double longitude = radians(LongitudeDegrees);

// solar panels
const double PanelBezel = 0.01; // m
const double PanelArea = (1.890 - PanelBezel) * (1.046 - PanelBezel);  // m^2
const double Efficiency = 0.204;
// https://www.enfsolar.com/pv/panel-datasheet/crystalline/55457

const int Group1Qty = 14;
const double Group1Azimuth = 0; // South

const int Group2Qty = 2;
const double Group2Azimuth = PI/2; // West


// loads
const double AirCon_kW = 5.6; // TODO: 2 stages?
const double Idle_kW = 0.55;   // reduce -- depends on PCs, fan, ...
const double HouseFan_kW = 0.25;  // depends on speed

// batteries
const double TMY_kWh = 75;
const int ChargerVolts = 240;
const double CMax_kWh = 7.2;
const double Powerwall_kWh = 13.5;
const int PowerwallMin_kWh = (20 - 4) * Powerwall_kWh / 100; // TODO: goes lower 
const double RTloss = 0.10; // depends on temperature, pWall fans

// References:
// https://www.scribd.com/document/725924868/7-1-Solar-Radiation-on-Inclined-Surfaces
// https://hal.science/hal-02175988/document
// https://emf-creaf.github.io/meteolandbook/solarradiation.html
// https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-on-a-tilted-surface

// TODO: 
//   model tree(s)
//   get outside temperature, wind to derate panel output
//   read Darkness for cloud cover adjust
//   panel temperature / wind model
// 
//   Powerwall level (API):
//      https://github.com/jasonacox/Powerwall-Dashboard/discussions/392
//      https://github.com/jasonacox/Powerwall-Dashboard/issues/387#issuecomment-1823104144
//      https://www.tesla.com/support/energy/powerwall/own/connecting-network
//      https://developer.tesla.com/docs/fleet-api/getting-started/what-is-fleet-api
//      https://github.com/vloschiavo/powerwall2
//      https://teslapWall_rkpfnt/   needs PW for AP gateway
//   TMY charge level
// 
//   PWall round trip adjust; PWall fans
//   air con - predict usage based on temperature forecast
//   adjust for house fan, recirc pump, ...
//   control Cmax outlet?


#define max(a,b) ((a) > (b) ? (a) : (b))

struct tm tm; 
#define day_of_year (tm.tm_yday)
#define orbit(day) (2 * PI * (day_of_year - (day)) / 365) 

double AST() { // Apparent Solar Time in hours
	time_t utc = time(NULL);
  gmtime_s(&tm, &utc);
	// double orb = 2 * PI * (day_of_year + 1 - 81) / 364;
	// double EoT = 9.87 * sin(2 * orb) - 7.53 * cos(orb) + 1.5 * sin(orb); // minutes

	double orb = 0.01720197 * (365.25 * (tm.tm_year - 100) + day_of_year) + 6.24004077;
	double EoT = -7.659 * sin(orb) + 9.863 * sin(2 * orb + 3.5932); // minutes
	return fmod((utc % (24 * 60 * 60)) / 3600. + longitude / radians(360 / 24) + EoT / 60 + 24, 24) - 12;  // -/+ 12 hours
}

double insolation(double solar_hour, double PanelAzimuth, double beta = asin(3./12)) { // solar hour +/-12; PanelAzimuth +West; beta = panel tilt
	// follows https://www.scribd.com/document/725924868/7-1-Solar-Radiation-on-Inclined-Surfaces

	double delta = radians(23.45) * sin(orbit(81));  // sun declination
	const double L = latitude;

	double h = 2 * PI * solar_hour / 24;  // hour angle	
	double alpha = asin(cos(L) * cos(delta) * cos(h) + sin(L) * sin(delta)); // solar altitude / elevation angle
  if (alpha < 0.01) return 0; // 0 when sun is below horizon

	const double A = 0.08;   // Altitude above sea level in km
	bool summer = day_of_year >= 80 && day_of_year <= 266;  // better smooth fn
	double r0 = summer ? 0.97 : 1.03;
	double r1 = summer ? 0.99 : 1.01;
	double rk = summer ? 1.02 : 1.00;
	double a0 = r0 * (0.4237 - 0.00821 * pow(6   - A, 2));
	double a1 = r1 * (0.5055 + 0.00595 * pow(6.5 - A, 2));
	double k =  rk * (0.2711 + 0.01858 * pow(2.5 - A, 2));
  double z = PI/2 - alpha; // zenith angle	
	double taub = a0 + a1 * exp(-k / cos(z)); // transmittance

	const double Isc = 1362; // solar "constant" W/m^2; varies +/- 5 with 11 year solar cycle
	double Ion = Isc * (1 + 0.033 * cos(orbit(0)));

  double phi = acos((sin(alpha) * sin(L) - sin(delta)) / cos(alpha) * cos(L)) * (h >= 0 ? 1 : -1); // solar azimuth +West
	double theta = acos(cos(alpha) * cos(fabs(phi - PanelAzimuth)) * sin(beta) + sin(alpha) * cos(beta)); // angle of incidence
	double Ib = Ion * taub * cos(theta);

	double taud = 0.2710 - 0.2939 * taub;
	double Ids = Ion * sin(alpha) * taud * (1 + cos(beta) / 2);

	const double rho = 0.1; // ground reflectance
	double Idg = (Ion * sin(alpha) * (taud + taub)) * rho * (1 - cos(beta)) / 2;

	return Ib + Ids + Idg;
}

double kW(double hour, double cellTemp = 45.7) {
	double tempDerate = (cellTemp - 25) * 0.0125; // temperature dependendence
	// https://www.enfsolar.com/pv/panel-datasheet/crystalline/55457

  double ic = Group1Qty * insolation(hour, Group1Azimuth) + Group2Qty * insolation(hour, Group2Azimuth);  // sum panel groups
	return ic * PanelArea * Efficiency * (1 - tempDerate) / 1000;
}

char response[256];

double getVal(const char* tag) {
  char* tagPos = strstr(response, tag);
  if (tagPos) {
    char* valPos = tagPos + strlen(tag) + 2; // ":
    return atof(valPos);
  }
  return 0;
}

bool airConOn() {
#ifdef ThermostatStatusURL
	FILE* resp = _popen("curl " ThermostatStatusURL " 2> NUL", "r");
  size_t got = fread(response, 1, sizeof response, resp);
  _pclose(resp);
	// {"temp":74.50,"tmode":2,"fmode":0,"override":0,"hold":0,"t_cool":72.00,"tstate":2,"fstate":1,"time":{"day":4,"hour":12,"minute":42},"t_type_post":0}
 	return  getVal("tstate") == 2;
#else
	return false;
#endif
}

void setAirConTargetTemp(int degF) {
#ifdef SetAirConTempURL
	char setTargetTempCmd[256];
	sprintf_s(setTargetTempCmd, sizeof setTargetTempCmd, "curl -d " SetAirConTempURL " > NUL", degF);
	system(setTargetTempCmd);  // TODO: returns success but no effect; ?needs auth??
#endif
}

int main() {
	const int pWallPercent = 27;  // TODO: enter to start
	bool isdst = true;  // TODO - set
	
	double pWall_kWh = pWallPercent * Powerwall_kWh / 100;

	while (1) {
	  printf("\033[2J\033[H"); // clear screen, home

	 double pWallPercent = 100 * pWall_kWh / Powerwall_kWh;

		double hour = AST();
		double kWnow = kW(hour);
		printf("%4.1f kW @ solar hour %.1f\n", kWnow, hour); 
		printf("\n kWh\n");

		double rest_of_day = 0;
		for (double hr = hour; hr <= 12; hr += 0.25) {
			double kw = kW(hr);  // TODO: temperature = f(t, Lo, Hi)
			if (hr > 0 && kw <= 0) break;		
			rest_of_day += kw; 
		}
		rest_of_day /= 4;
		if (rest_of_day > 0) 
			printf("%4.1f solar rest of day\n", rest_of_day);

	  double sundownHr;
		for (double hr = 4; hr <= 12; hr += 0.25) 
			if (kW(hr) <= Idle_kW) {
				sundownHr = hr;
				break;
			}

		double idleToSundown = hour < sundownHr ? (sundownHr - hour) * Idle_kW : 0;
		double overnight_kWh = (12 - sundownHr) * Idle_kW;  // midnite to morn

		double idleToMidnite = max(0, (12 - hour - isdst)) * Idle_kW; // TODO: accurate midnite
    double pWallToMidnite = idleToMidnite + PowerwallMin_kWh;
	  double excessToMidnite = pWall_kWh + rest_of_day - pWallToMidnite;
	  if (rest_of_day > 0 && pWallToMidnite > pWall_kWh) 
			printf("%4.1f PW %.0f%% to %.0f%% til midnite\n", -(pWallToMidnite - pWall_kWh), pWallPercent, 100 * pWallToMidnite / Powerwall_kWh);		
	
		double pWallToMorn = pWallToMidnite + overnight_kWh;
		double excessToMorn = pWall_kWh + rest_of_day - pWallToMorn;
		if (rest_of_day > 0 && pWallToMorn > pWall_kWh) 
			printf("%4.1f PW %.0f%% to %.0f%% for overnight\n", -(pWallToMorn - pWallToMidnite), 100 * pWallToMidnite / Powerwall_kWh, 100 * pWallToMorn / Powerwall_kWh);

		if (excessToMorn > 0) {
			printf("%4.1f excess to morn\n", excessToMorn);
			printf("\n%3d/%d%% TMY charge OR\n", (int)(100 * excessToMorn / TMY_kWh), (int)(100 * excessToMidnite / TMY_kWh));
			printf("%.1f/%.1f hours air con\n", excessToMorn / AirCon_kW, excessToMidnite / AirCon_kW);
		} else {
			if (excessToMidnite < 0) 
				printf("%4.1f from grid < midnite\n", -excessToMidnite);
			printf("%4.1f from grid > midnite\n", -excessToMorn - excessToMidnite);
			setAirConTargetTemp(excessToMidnite < 0 ? 80 : 78);
		}

		bool airCon = airConOn();
		if (airCon)
			printf("%4.1f air con %.1f --> %.0f deg\n", -AirCon_kW, getVal("temp"), getVal("t_cool"));
		else if ((kWnow - Idle_kW) > 5 * ChargerVolts / 1000. && pWallPercent < 25)
		  printf("%4.0f Amp no grid TMY charge\n", (kWnow - Idle_kW) * 1000 / ChargerVolts); // to avoid using grid
		double homeKw = Idle_kW + (airCon ? AirCon_kW : 0);  // TODO: plus TMY, CMax, house fan, ...
  	static double prevHour = hour;			
		pWall_kWh += (hour - prevHour) * (kWnow - homeKw);
		prevHour = hour;
		if (pWall_kWh > Powerwall_kWh) pWall_kWh = Powerwall_kWh;
		if (pWall_kWh < PowerwallMin_kWh) pWall_kWh = PowerwallMin_kWh; // TODO -- can be less

		Sleep(15 * 60 * 1000);
	} 

  return 0;
}


#if 0 // alternate calculations -- some with more accurate constants??

double Ea() { // = Io; at top of atmosphere
  const double Esc = 1362; // solar constant W / m^2
	double b = orbit(0);
  return Esc * 1.00011 + 0.034221 * cos(b) + 0.00128/sin(b) + 0.00719 * cos(2 * b) + 0.00007 * sin(2 * b);
}

// solar_hour: from solar noon
double insolation_bad(double solar_hour, double PanelAzimuth, double Z = asin(3./12), double reflectance = 0.1) {  // Z = collector tilt
	double H = 2 * PI * solar_hour / 24;  // hour angle
  double theta = radians(23.45) * sin(orbit(81));  // Declination of the sun, 0 at spring and fall equinoxes
  double beta = asin(cos(latitude) * cos(theta) * cos(H) + sin(latitude) * sin(theta));  // Sun Altitude angle
	  
  double A = 1160 + 75 * sin(orbit(275)); // Apparent extraterrestrial flux
	double sun_azimuth = A * sin(cos(theta) * sin(H) / cos(beta));
  if (cos(H) < tan(theta) / tan(latitude)) // never?
		sun_azimuth = PI - sun_azimuth;

  double k = 0.174 + 0.035 * sin(orbit(100));  // Optical depth
	double m = sqrt(pow(708 * sin(beta), 2) + 1417) - 708 * sin(beta); // = fabs(1 / sin(beta));    // air mass ratio
  double Ib = A * exp(-k * m);
  double Ibc = Ib * (cos(sun_azimuth - PanelAzimuth) * cos(beta) * sin(Z) + sin(beta) * cos(Z));  // direct beam radiation
 
	double C = 0.095 + 0.04 * sin(orbit(100));  // sky diffuse factor
  double Idc = Ib * C * (1 + cos(Z)) / 2;  // diffuse radiation
  double Irc = reflectance * Ib * (C + sin(beta))*(1 - cos(Z)) / 2;  // reflected radiation
  
	double Ic = Ibc + Idc + Irc; // Total insolation
  return Ic;  // W / m^2
}

    // at solar noon
  double Sdecl = 23.45 * sin(orbit(81));
  double Salt = 90

  // Insolation
  double A = 1160 + 75 * sin(orbit(275)); // Apparent extraterrestrial flux
  double k = 0.174 + 0.035 * sin(orbit(100)); // Optical depth
  double C = 0.095 + 0.04 * sin(orbit(100));

double getsunpos(LAT, LON, year, month, day, hour, &azimuth, &altitude, &sunstate) {
	
	  // julian date
	d = 367*year - floor((7*(year + floor((month+9)/12)))/4)
				 + floor((275*month)/9) + day - 730530;
		
	w = 4.9382415669097640822661983551248 
			 + .00000082193663128794959930855831205818* d; // (longitude of perihelion)
	  a = 1.000000                           ;//    (mean distance, a.u.)
	  e = 0.016709 - .000000001151        * d ;//   (eccentricity)
	  M = 6.2141924418482506287494932704807 
	  	 + 0.017201969619332228715501852561964 * d ;//   (mean anomaly)
		
		
	oblecl = 0.40909295936270689252387465029835 
						- .0000000062186081248557962825791102081249 * d  ;// obliquity of the ecliptic
		
	latitude = w + M; // sun's mean longitude
		
	E = M + e * sin(M) * (1 + e * cos(M));
		
	x = cos(E) - e;
        y = sin(E) * sqrt(1 - e * e);
	    
	  r = sqrt(x*x + y*y);
	  v = atan2( y, x )  ;
	    
	  lon = v + w;
	  
	  x = r * cos(lon);
	  y = r * sin(lon);
	  z = 0.0;
	  
	  xequat = x;
	  yequat = y * cos(oblecl) + z * sin(oblecl);
	  zequat = y * sin(oblecl) + z * cos(oblecl);
	
		RA   = atan2( yequat, xequat );
	  Decl = asin( zequat / r );
	
		RAhours = r2d(RA)/15;
		
	  GMST0 = r2d( latitude + pi() ) / 15;//
	  SIDTIME = GMST0 + hour + rad2deg(LON)/15;
	  
		HA = deg2rad((SIDTIME - RAhours) *15);
		
		x = cos(HA) * cos(Decl);
	  y = sin(HA) * cos(Decl);
	  z = sin(Decl);
	  
	  xhor = x * cos(pi()/2 - LAT) - z * sin(pi()/2 - LAT);
	  yhor = y;
	  zhor = x * sin(pi()/2 - LAT) + z * cos(pi()/2 - LAT);
	  
	  azimuth = rad2deg(atan2( yhor, xhor ) + pi());
	  altitude = rad2deg(asin( zhor ));
	  	  
	}	

	LAT = deg2rad(LatitudeDegrees);
	LON = deg2rad(LongitudeDegrees);
	
	year = gmdate("Y");
	month = gmdate("m");
	day = gmdate("d");
	hour = gmdate("H") + (gmdate("i") / 60);

  // get current position
	getsunpos(LAT, LON, year, month, day, hour, azimuth, altitude, sunstate);

	modtilt = 40; //Solar panel's tilt angle
	modazi = 200;//Solar panel's azimuth
	modsufrace=1.88512; //Solar panel's surface in sq. meters

 	//The intensity of the direct component of sunlight throughout each day can be determined as a function of air mass 
	//http://pveducation.org/pvcdrom/properties-of-sunlight/air-mass#formula      
	airmass = 1/cos((90-altitude) * 4 * asin(1) /360); 

  //Sincident is the intensity on a plane perpendicular to the sun's rays in units of kW/m2 and AM is the air mass. The value of 1.353 kW/m2 is the solar constant and the number 0.7 arises from the fact that about 70% of the radiation incident on the atmosphere is transmitted to the Earth. The extra power term of 0.678 is an empirical fit to the observed data and takes into account the non-uniformities in the atmospheric layers.
	Sincident=(1.353*pow(0.7,pow(airmass,0.678)));

	//A module that directly faces the sun so that the incoming rays are perpendicular to the module surface has the module tilt equal to the sun's zenith angle (90 - a = ß), and the module azimuth angle equal to the sun's azimuth angle 
	//Solar tubes are an example where the module azimuth can be treated as hooked to the solar's azimuth foe a -90/+90 degree angle
	//Comment out the following if you are on FLAT panel, not solar tubes:
	if(azimuth>(modazi-80) && azimuth<(modazi+80)) {
		modazi=azimuth;
	}

  fraction = cos(altitude*4*asin(1)/360)*sin(modtilt*4*asin(1)/360)*cos(azimuth*4*asin(1)/360-modazi*4*asin(1)/360)+sin(altitude*4*asin(1)/360)*cos(modtilt*4*asin(1)/360);

	// kW/m² light intensity on the module * module's surface
	Smodule = Sincident * fraction * modsufrace * 1000;

	if (Smodule<0) Smodule=0;

	if ( altitude < 0 ) {
		altitude=0;
   	azimuth=0;
		Smodule=0;
	}

	eqtime(jd) {
	if (nargs() < 1 ) {cat("USAGE: eqtime(jd)\n"); return()}
	jdc = (jd - 2451545.0)/36525.0
	sec = 21.448 - jdc*(46.8150 + jdc*(0.00059 - jdc*(0.001813)))
	e0 = 23.0 + (26.0 + (sec/60.0))/60.0 
	ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
	oblcorr = e0 + 0.00256 * cos(radians(125.04 - 1934.136 * jdc)) 
	y = (tan(radians(oblcorr)/2))^2
	l0 = 280.46646 + jdc * (36000.76983 + jdc*(0.0003032))
	l0 = (l0-360*(l0%/%360))%%360
	rl0 = radians(l0)
	gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
	gmas = radians(gmas)
	EqTime = y*sin(2*rl0)-2.0*ecc*sin(gmas)+4.0*ecc*y*sin(gmas)*cos(2*rl0)-
		0.5*y^2*sin(4*rl0)-1.25*ecc^2*sin(2*gmas)
	return(degrees(EqTime)*4)
}

declination(jd) {
    if (nargs() < 1) {
        cat("USAGE: declination(jd) \n jd = Julian day \n")
        return()
    }
  ;  // Julian Centuries (Meeus, Astronomical Algorithms 1999. (24.1))
    T = (jd - 2451545)/36525.0
  ;  // mean obliquity of the ecliptic (21.2)
    epsilon = (23+26/60.0+21.448/3600.0) - (46.8150/3600.0)*T - 
    			(0.00059/3600.0)*T^2 + (0.001813/3600.0)*T^3
  ;  // geometric mean longitude of the sun (24.2)
    L0 = 280.46645 + 36000.76983*T + 0.0003032*T^2
  ;  // L0 = (L0 - 360 * (L0%/%360))%%360
  ;  // mean anomaly of the Sun (24.3)
    M = 357.52910 + 35999.05030*T - 0.0001559*T^2 - 0.00000048*T^3
  ;  // eccentricity of the Earth's orbit (24.4)
    e = 0.016708617 - 0.000042037*T - 0.0000001236*T^2
  ;  // Sun's equation of center
    C = (1.914600 - 0.004817*T - 0.000014*T^2)*sin(radians(M)) + 
        (0.019993 - 0.000101*T)*sin(2*radians(M)) +
        0.000290*sin(3*radians(M))
  ;  // Sun's true longitude
    Theta = L0 + C
  ;  // Sun's true anomaly
    v = M + C
  ;  // Sun's Radius Vector (24.5)
  ;  // R = (1.000001018*(1-e^2))/(1 + e*cos(radians(v)))
  ;  //  Longitude of the ascending node of the moon
    Omega = 125.04452 - 1934.136261*T +0.0020708*T^2 +(T^3)/450000
  ;  // Apparent longitude of the sun
    lambda = Theta - 0.00569 - 0.00478*sin(radians(Omega))
  ;  // Sun's declination (24.7)
    delta = asin(sin(radians(epsilon)) * sin(radians(lambda)))
    return(degrees(delta))
}

#endif
