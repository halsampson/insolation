// insolation.cpp

// Solar panel production and Powerwall charge calculations

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
const double A = 0.08;   // Altitude above sea level in km

// solar panels
const int Group1Qty = 14;
const double Group1Azimuth = 0; // South

const int Group2Qty = 2;
const double Group2Azimuth = PI/2; // West

// typical panel specs
const double PanelBezel = 0.01; // m
const double PanelArea = (1.890 - PanelBezel) * (1.046 - PanelBezel);  // m^2
const double Efficiency = 0.204;
const double TempDerate = 1.25; // % / degC
// https://www.enfsolar.com/pv/panel-datasheet/crystalline/55457

// batteries
const double TMY_kWh = 75;
const int ChargerVolts = 240;
const double CMax_kWh = 7.2;

const double PWallMax_kWh = 13.5;  // Powerwall 3
const double PWallMin_kWh = (20 - 1) * PWallMax_kWh / 100; // TODO: why much lower than reserve early morning?
const double RTloss = 0.10; // depends on temperature, pWall fans  TODO: account for RTloss

// loads
const double AirCon_kW = 4.25;  // TODO: 2 stages?
const double Idle_kW = 0.6;   // reduce -- fridge depends on interior temp; PCs, ...
const double HouseFan_kW = 0.25;  // depends on speed

// weather
#define degC(F)  (((F) - 32) * 5 / 9.)
const double ForecastHighTemp = degC(88);
const double ForecastLowTemp  = degC(62);
const double WindSpeed = 8; // mph when hottest
// TODO: get ambient, wind each hour from weather forecast; also cloud cover %


// References:
// https://www.scribd.com/document/725924868/7-1-Solar-Radiation-on-Inclined-Surfaces
// https://hal.science/hal-02175988/document
// https://emf-creaf.github.io/meteolandbook/solarradiation.html
// https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-on-a-tilted-surface

// TODO: 
//   model tree(s)
//   get outside temperature, wind to derate panel output
//     better panel temperature / wind cooling model
//   read Darkness for cloud cover adjust
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

time_t utc;
struct tm tm; 
#define dayOfYear tm.tm_yday
#define orbit(day) (2 * PI * (dayOfYear - (day)) / 365) 
#define solarHr(hr) (fmod(hr + 24 + 12, 24) - 12)
double EoT;

const double L = latitude;
double delta; // sun declination
double a0, a1, k;  // air mass coefficients
double Ion;

void setDay() {
	utc = time(NULL);
  gmtime_s(&tm, &utc);

	static int prevDayOfYear = -1;
	if (dayOfYear == prevDayOfYear) return;

	// double orb = 2 * PI * (dayOfYear + 1 - 81) / 364;
	// double EoT = 9.87 * sin(2 * orb) - 7.53 * cos(orb) + 1.5 * sin(orb); // minutes

	double orb = 0.01720197 * (365.25 * (tm.tm_year - 100) + dayOfYear) + 6.24004077;
	double EoT = -7.659 * sin(orb) + 9.863 * sin(2 * orb + 3.5932); // Equation of Time minutes

	delta = radians(23.45) * sin(orbit(81));  
  bool summer = dayOfYear >= 80 && dayOfYear <= 266;  // better smooth fn
	double r0 = summer ? 0.97 : 1.03;
	double r1 = summer ? 0.99 : 1.01;
	double rk = summer ? 1.02 : 1.00;

	a0 = r0 * (0.4237 - 0.00821 * pow(6   - A, 2));
	a1 = r1 * (0.5055 + 0.00595 * pow(6.5 - A, 2));
	k =  rk * (0.2711 + 0.01858 * pow(2.5 - A, 2));

	const double Isc = 1362; // solar "constant" W/m^2; varies +/- 5 with 11 year solar cycle
	Ion = Isc * (1 + 0.033 * cos(orbit(0)));
}

double apparentSolarTime() { // Apparent Solar Time in hours
	return solarHr((utc % (24 * 60 * 60)) / 3600. + longitude / radians(360 / 24) + EoT / 60 - 12);  // -/+ 12 hours from solar noon
}

double insolation(double solar_hour, double PanelAzimuth, double beta = asin(3./12)) { // solar hour +/-12; PanelAzimuth +West; beta = panel tilt
	// follows https://www.scribd.com/document/725924868/7-1-Solar-Radiation-on-Inclined-Surfaces
	double h = 2 * PI * solar_hour / 24;  // hour angle	
	double alpha = asin(cos(L) * cos(delta) * cos(h) + sin(L) * sin(delta)); // solar altitude / elevation angle
  if (alpha < 0.01) return 0; // 0 when sun is below horizon

  double z = PI/2 - alpha; // zenith angle	
	double taub = a0 + a1 * exp(-k / cos(z)); // transmittance

  double phi = acos((sin(alpha) * sin(L) - sin(delta)) / cos(alpha) * cos(L)) * (h >= 0 ? 1 : -1); // solar azimuth +West
	double theta = acos(cos(alpha) * cos(fabs(phi - PanelAzimuth)) * sin(beta) + sin(alpha) * cos(beta)); // angle of incidence
	double Ib = Ion * taub * cos(theta);

	double taud = 0.2710 - 0.2939 * taub;
	double Ids = Ion * sin(alpha) * taud * (1 + cos(beta) / 2);

	const double rho = 0.2; // ground reflectance
	double Idg = (Ion * sin(alpha) * (taud + taub)) * rho * (1 - cos(beta)) / 2;

	return Ib + Ids + Idg;
}

double kW(double hour) { 
	double mainInsolation = insolation(hour, Group1Azimuth);
  double ic = Group1Qty * mainInsolation  + Group2Qty * insolation(hour, Group2Azimuth);  // sum panel groups

	double ambientTemp = ForecastHighTemp * fabs(solarHr(hour + 12 - 3)) / 12 
		                 + ForecastLowTemp  * fabs(solarHr(hour - 3)) / 12; // high at solar hour 3
	double tempRise = (25.7 * mainInsolation / 800) * max(0, 1 - WindSpeed / 30);  // TODO: better estimate: f(wind direction), radiative 
  double tempDerate = (ambientTemp + tempRise - 25) * TempDerate / 100 ; // temperature dependence
	// https://www.enfsolar.com/pv/panel-datasheet/crystalline/55457

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
	system(setTargetTempCmd);  // TODO: returns success but no effect; ?unlock?
#endif
}

int main() {
	while (1) {
		printf("\nPowerwall level   %%\b\b\b");
		int pWallPercent = 10 * (_getche() - '0') + _getche() - '0'; 
		if (pWallPercent > 100) pWallPercent = 100; // letters -> 100%
		double pWall_kWh = pWallPercent * PWallMax_kWh / 100;

		bool isdst = true;  // TODO - set
	
		while (1) {
			printf("\033[2J\033[H"); // clear screen, home

			double pWallPercent = 100 * pWall_kWh / PWallMax_kWh;

			setDay();
			double hour = apparentSolarTime();
			double kWnow = kW(hour);
			printf("%4.1f kW @ solar hour %.1f\n", kWnow, hour); 
			printf("\n kWh\n");

			double generated = 0;
			for (double hr = -8; hr <= hour; hr += 0.25) 
				generated += kW(hr); 
			generated /= 4;
			if (generated > 0) 
				printf("%4.1f generated today\n", generated); // TODO: tree, clouds adjust

			double rest_of_day = 0;
			for (double hr = hour; hr < 8; hr += 0.25) {
				double kw = kW(hr);
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

			double idleToSundown = max(0, sundownHr - hour) * Idle_kW;	
			printf("%4.1f idle to sundown\n", -idleToSundown);

			double idleSundownToMidnite = max(0, (12 - isdst - max(hour, sundownHr))) * Idle_kW; // TODO: accurate midnite
			double pWallSundownToMidnite = idleSundownToMidnite + PWallMin_kWh;
			double excessToMidnite = pWall_kWh + rest_of_day - idleToSundown- pWallSundownToMidnite;
			if (excessToMidnite >= 0 && pWall_kWh < pWallSundownToMidnite) 
				printf("%4.1f PW %.0f%% to %.0f%% sundown to midnite\n", pWall_kWh - pWallSundownToMidnite, pWallPercent, 100 * pWallSundownToMidnite / PWallMax_kWh);		

			double overnight_kWh = (12 + isdst - sundownHr) * Idle_kW;  // midnite to morn sunup
			double pWallToMorn = pWallSundownToMidnite + overnight_kWh;
			double excessToMorn = excessToMidnite - overnight_kWh;
			if (excessToMorn >= 0 && pWall_kWh < pWallToMorn)
				printf("%4.1f PW %.0f%% to %.0f%% for overnight\n", -overnight_kWh, 100 * max(pWall_kWh, pWallSundownToMidnite) / PWallMax_kWh, 100 * pWallToMorn / PWallMax_kWh);

			if (excessToMorn > 0) {
				printf("%4.1f excess to next morning\n", excessToMorn);
				printf("\n%.1f/%.1f hours air con OR\n", excessToMorn / AirCon_kW, excessToMidnite / AirCon_kW);
				printf("%3d/%d%% TMY charge\n", (int)(100 * excessToMorn / TMY_kWh), (int)(100 * excessToMidnite / TMY_kWh));
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
			if (hour > prevHour) // avoid midnite
				pWall_kWh += (hour - prevHour) * (kWnow - homeKw);
			prevHour = hour;
			if (pWall_kWh > PWallMax_kWh) pWall_kWh = PWallMax_kWh;
			if (pWall_kWh < PWallMin_kWh) pWall_kWh = PWallMin_kWh; 

			int sleep = 15 * 60 * 1000 / 100;
			while (!_kbhit() && sleep--)
				Sleep(100);
			if (_kbhit()) break;
		} 
	}

  return 0;
}

