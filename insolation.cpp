// insolation.cpp

// Solar panel production and Powerwall charge calculations

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <windows.h>

#include "private.txt"  // #define LatitudeDegrees and LongitudeDegrees

const double PI = 3.14159265359;
#define radians(degrees)  (PI * (degrees) / 180)

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
const int PowerwallMin_kWh = (20 - 14) * Powerwall_kWh / 100; // TODO: why lower  than reserve??
const double RTloss = 0.10; // depends on temperature, pWall fans  TODO: account for RTloss

const double latitude = radians(LatitudeDegrees);
const double longitude = radians(LongitudeDegrees);

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

struct tm tm; 
#define day_of_year tm.tm_yday
#define orbit(day) (2 * PI * (day_of_year - (day)) / 365) 

double AST() { // Apparent Solar Time in hours
	time_t utc = time(NULL);
  gmtime_s(&tm, &utc);
	// double orb = 2 * PI * (day_of_year + 1 - 81) / 364;
	// double EoT = 9.87 * sin(2 * orb) - 7.53 * cos(orb) + 1.5 * sin(orb); // minutes

	double orb = 0.01720197 * (365.25 * (tm.tm_year - 100) + day_of_year) + 6.24004077;
	double EoT = -7.659 * sin(orb) + 9.863 * sin(2 * orb + 3.5932); // Equation of Time minutes
	return fmod((utc % (24 * 60 * 60)) / 3600. + longitude / radians(360 / 24) + EoT / 60 + 24, 24) - 12;  // -/+ 12 hours from solar noon
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
	printf("Powerwall level   %%\b\b\b");
  int pWallPercent = 10 * (_getche() - '0') + _getche() - '0'; 
	if (pWallPercent > 100) pWallPercent = 100; // letters -> 100%
	double pWall_kWh = pWallPercent * Powerwall_kWh / 100;

	bool isdst = true;  // TODO - set
	
	while (1) {
	  printf("\033[2J\033[H"); // clear screen, home

	  double pWallPercent = 100 * pWall_kWh / Powerwall_kWh;

		double hour = AST();
		double kWnow = kW(hour);
		printf("%4.1f kW @ solar hour %.1f\n", kWnow, hour); 
		printf("\n kWh\n");

		double rest_of_day = 0;
		for (double hr = hour; hr < 12; hr += 0.25) {
			double kw = kW(hr);  // TODO: temperature = f(t, Lo, Hi)
			if (hr > 0 && kw <= 0) break;		
			rest_of_day += kw; 
		}
		rest_of_day /= 4;
		rest_of_day *= 1 - 0.08;  // TODO: panel temperature, tree adjust - compare to actual
		if (rest_of_day > 0) 
			printf("%4.1f solar rest of day\n", rest_of_day);

	  double sundownHr;
		for (double hr = 4; hr <= 12; hr += 0.25) 
			if (kW(hr) <= Idle_kW) {
				sundownHr = hr;
				break;
			}

		double idleToSundown = hour < sundownHr ? (sundownHr - hour) * Idle_kW : 0;		
		double idleToMidnite = (12 - isdst - hour - sundownHr) * Idle_kW; // TODO: accurate midnite
		double overnight_kWh = (12 + isdst - sundownHr) * Idle_kW;  // midnite to morn sunup

    double pWallToMidnite = idleToMidnite + PowerwallMin_kWh;
	  double excessToMidnite = pWall_kWh + rest_of_day - pWallToMidnite;
	  if (rest_of_day > 0 && pWallToMidnite > pWall_kWh) 
			printf("%4.1f PW %.0f%% to %.0f%% til midnite\n", -(pWallToMidnite - pWall_kWh), pWallPercent, 100 * pWallToMidnite / Powerwall_kWh);		
	
		double pWallToMorn = pWallToMidnite + overnight_kWh;
		double excessToMorn = pWall_kWh + rest_of_day - pWallToMorn;
		if (rest_of_day > 0 && pWallToMorn > pWall_kWh) 
			printf("%4.1f PW %.0f%% to %.0f%% for overnight\n", -(pWallToMorn - pWallToMidnite), 100 * max(pWall_kWh, pWallToMidnite) / Powerwall_kWh, 100 * pWallToMorn / Powerwall_kWh);

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
		if (hour > prevHour) // avoid midnite
		  pWall_kWh += (hour - prevHour) * (kWnow - homeKw);
		prevHour = hour;
		if (pWall_kWh > Powerwall_kWh) pWall_kWh = Powerwall_kWh;
		if (pWall_kWh < PowerwallMin_kWh) pWall_kWh = PowerwallMin_kWh; 
		Sleep(15 * 60 * 1000);
	} 

  return 0;
}

