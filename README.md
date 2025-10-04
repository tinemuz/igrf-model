# IGRF Model Library

A lightweight Java library for calculating Earth's magnetic field using the International Geomagnetic Reference Field (IGRF) model. Essential for navigation applications requiring accurate magnetic declination, inclination, and field intensity calculations.

## Features

- **Ready to use** - IGRF coefficients included, no configuration needed
- **Zero allocation** after initialization - ideal for real-time systems
- **Thread-safe** with per-thread workspace caching
- **High accuracy** using IGRF spherical harmonic model (degree 13)
- **No external dependencies** - pure Java implementation
- **Small footprint** - single class with embedded coefficients
- **WGS-84 support** for geodetic coordinates

## Installation

### Maven
```xml
<dependency>
    <groupId>com.github.tinemuz</groupId>
    <artifactId>igrf-model</artifactId>
    <version>1.0-SNAPSHOT</version>
</dependency>
```

### Gradle
```gradle
implementation 'com.github.tinemuz:igrf-model:1.0-SNAPSHOT'
```

## Usage

### Quick Start

The library provides two primary methods: `declinationDeg()` for simple magnetic variation calculations and `compute()` for full field data.

**Initialize at startup:**
```java
IGRFModel.preload(); // Recommended to call at application startup
```

### Calculate Magnetic Declination

Use this for converting between true and magnetic headings:

```java
double declination = IGRFModel.declinationDeg(
    40.7128,                    // Latitude (degrees, north positive)
    -74.0060,                   // Longitude (degrees, east positive)
    100.0,                      // Altitude (meters above MSL)
    System.currentTimeMillis()  // Time (epoch milliseconds UTC)
);

// Convert headings
double magneticHeading = trueHeading - declination;
double trueHeading = magneticHeading + declination;
```

### Get Full Magnetic Field

For all field components and derived quantities:

```java
Result field = IGRFModel.compute(latitude, longitude, altitude, epochMillis);

// Access components (all in nanoteslas, except angles in degrees)
field.xNorthNt       // North component
field.yEastNt        // East component  
field.zDownNt        // Down component
field.hHorizontalNt  // Horizontal intensity
field.fTotalNt       // Total field intensity
field.declinationDeg // Declination (magnetic variation)
field.inclinationDeg // Inclination (dip angle)
```

### Common Applications

**Aviation - Flight Planning:**
```java
double variation = IGRFModel.declinationDeg(lat, lon, cruiseAlt, departureTime);
double magneticCourse = trueCourse - variation;
```

**Marine - Compass Navigation:**
```java
Result field = IGRFModel.compute(vesselLat, vesselLon, 0, System.currentTimeMillis());
double trueHeading = compassHeading + field.declinationDeg + compassDeviation;
```

**Survey - Magnetic Anomaly Detection:**
```java
Result field = IGRFModel.compute(surveyLat, surveyLon, surveyAlt, surveyTime);
double anomaly = measuredField - field.fTotalNt;
```

**Drones - Magnetometer Calibration:**
```java
Result field = IGRFModel.compute(droneLat, droneLon, droneAlt, System.currentTimeMillis());
calibrator.setReference(field.xNorthNt, field.yEastNt, field.zDownNt);
```

## API Reference

### Methods

**`IGRFModel.preload()`**
- Pre-loads IGRF coefficients from classpath
- Recommended to call at application startup
- Validates model and prints warning if time is >5 years beyond latest epoch

**`IGRFModel.declinationDeg(latitude, longitude, altitude, epochMillis)`**
- Returns magnetic declination in degrees (east positive, west negative)
- Fastest method when only declination is needed
- Range: [-180°, +180°)

**`IGRFModel.compute(latitude, longitude, altitude, epochMillis)`**
- Returns `Result` object with all magnetic field components
- Components in nanoteslas (nT), angles in degrees

### Result Fields

| Field | Unit | Description |
|-------|------|-------------|
| `xNorthNt` | nT | North component |
| `yEastNt` | nT | East component |
| `zDownNt` | nT | Down component |
| `hHorizontalNt` | nT | Horizontal intensity |
| `fTotalNt` | nT | Total field intensity |
| `declinationDeg` | degrees | Magnetic declination/variation |
| `inclinationDeg` | degrees | Magnetic inclination/dip |

## Coordinate Systems

**Input (WGS-84 Geodetic):**
- Latitude: -90° (South) to +90° (North)
- Longitude: -180° to +180° or 0° to 360°
- Altitude: meters above WGS-84 ellipsoid

**Output (NED Frame):**
- X: Geographic north
- Y: Geographic east
- Z: Downward (toward Earth's center)

All field components are in nanoteslas (nT).

## Technical Details

### IGRF Model

The library implements the International Geomagnetic Reference Field using spherical harmonic expansion to degree and order 13. The magnetic scalar potential V at position (r, θ, φ) is:

```
V = a Σ(a/r)^(n+1) Σ[g_n^m cos(mφ) + h_n^m sin(mφ)]P_n^m(cos θ)
```

Where:
- a = 6371.2 km (reference Earth radius)
- g_n^m, h_n^m are Gauss coefficients (Schmidt quasi-normalized)
- P_n^m are associated Legendre functions
- N = 13 (maximum degree and order)

The magnetic field components are derived as B = -∇V, with additional transformations from geocentric to geodetic coordinates.

### Key Algorithms

**Coordinate Conversion:** Geodetic (WGS-84) to geocentric coordinates using ellipsoid parameters.

**Time Interpolation:** Linear interpolation between IGRF epochs, with secular variation extrapolation up to 5 years beyond the latest epoch.

**Field Computation:** Stable recurrence relations for associated Legendre functions with Schmidt normalization.

**Optimizations:** Thread-local workspace caching, trigonometric series caching, and zero-allocation design for high-performance repeated computations.

## Model Accuracy & Typical Values

### Typical Field Values by Location

| Location | Total Field (nT) | Declination | Inclination |
|----------|------------------|-------------|-------------|
| Equator | 25,000-40,000 | ±20° | ~0° (horizontal) |
| Mid-latitudes | 45,000-55,000 | ±30° | 50-70° |
| Polar regions | 55,000-65,000 | Variable | ~90° (vertical) |

### Accuracy (IGRF-13, 2015-2025)

| Quantity | Typical Accuracy |
|----------|------------------|
| Declination | ±0.5° (±1° near poles) |
| Inclination | ±0.5° (±1° near equator) |
| Total Intensity | ±100-200 nT |
| Components (X, Y) | ±100-150 nT |
| Component (Z) | ±150-250 nT |

**Important:** Accuracy degrades beyond 5 years past the latest epoch. Local magnetic anomalies can cause deviations of thousands of nT. Update coefficients regularly for best results.

## Performance

**Initialization:**
- First call: ~50-100ms (loads coefficients)
- Memory: ~200KB coefficients + ~20KB per thread

**Runtime:**
- Single evaluation: ~5-15 microseconds
- Allocations: Zero per call after initialization
- Thread safety: Lock-free via ThreadLocal workspace

## Time Handling

The library converts epoch milliseconds (UTC) to decimal years internally:

```java
// Current time
long now = System.currentTimeMillis();

// Specific date
LocalDateTime dt = LocalDateTime.of(2025, 1, 15, 12, 0);
long epochMillis = dt.toInstant(ZoneOffset.UTC).toEpochMilli();
```

**Time periods:**
- Historical: Linear interpolation between IGRF epochs
- Future: Secular variation applied up to +5 years (then clamped with warning)

## IGRF Coefficients

The library comes with **embedded IGRF coefficients** included in the JAR. Users don't need to provide or manage coefficient files - everything works out of the box.

The current version includes IGRF-14 coefficients covering epochs 1900-2030. When new IGRF models are released, simply update the library version to get the latest coefficients.

## FAQ

**What is magnetic declination?**  
The angle between true north (geographic) and magnetic north. It varies by location and time.

**How do I convert between true and magnetic headings?**  
Magnetic heading = True heading - Declination  
True heading = Magnetic heading + Declination

**What about compass deviation?**  
This library provides magnetic declination (variation). Compass deviation from local interference must be determined separately through calibration.

**Thread safety?**  
Fully thread-safe with per-thread workspace caching via ThreadLocal.

**Historical data support?**  
Yes, the model covers 1900-2025 with linear interpolation between epochs.

## Troubleshooting

**"requested time is X years beyond latest epoch"**  
The calculation date is more than 5 years beyond the latest IGRF epoch included in the library. Results are clamped to +5 years with secular variation applied. Consider updating to a newer version of the library if available.

**"IGRF coefficients file not found"**  
This indicates a packaging issue with the JAR. The coefficients should be automatically included. Verify you're using an official release or rebuild the library properly.

**Performance issues**  
Call `preload()` at startup and reuse threads when possible to benefit from cached workspaces.

## References

- [IGRF Model Documentation](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
- [NOAA Geomagnetic Calculators](https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml)
- [WGS-84 Specification](https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84)

## License

[Specify your license here]

## Contributing

Contributions welcome! Please ensure code follows existing conventions, all tests pass, and documentation is updated.

---

**Version 1.0-SNAPSHOT** - IGRF-14 model support (epochs 1900-2030)
