# IGRF Model Library

A lightweight Java library for calculating Earth's magnetic field using the International Geomagnetic Reference Field (IGRF) model. The IGRF is a model of the Earth's main magnetic field that is updated every 5 years. This library provides accurate magnetic declination, inclination, and field intensity calculations essential for navigation applications.

The library includes embedded IGRF-14 coefficients (epochs 1900-2030) in the JAR. The implementation uses spherical harmonic expansion to degree and order 13. It has no external dependencies and supports WGS-84 geodetic coordinates.

No need to source and maintain your own coefficient files.

> Note that this library uses IGRF-14 coefficients (epochs 1900-2030)

## Install

### Maven
```xml
<dependency>
    <groupId>io.github.tinemuz</groupId>
    <artifactId>igrf-model</artifactId>
    <version>0.0.1-SNAPSHOT</version>
</dependency>
```

### Gradle
```gradle
implementation 'io.github.tinemuz:igrf-model:0.0.1-SNAPSHOT'
```

> **Note:** Snapshot versions are published automatically to Maven Central snapshots repository. For snapshots, add the repository:
> ```xml
> <repositories>
>     <repository>
>         <id>sonatype-snapshots</id>
>         <url>https://central.sonatype.com/repository/maven-snapshots/</url>
>         <snapshots>
>             <enabled>true</enabled>
>         </snapshots>
>     </repository>
> </repositories>
> ```

### Copy Into Your Project

If you prefer not to add a dependency, you can copy the source files directly into your project. The library consists of just two files:

- `IGRFModel.java` - The main computation class
- `igrfcoeffs.txt` - IGRF-14 coefficients (must be placed in your resources' folder)

## Example

It is recommended to preload the model at application startup:

```java
IGRFModel.preload(); // Recommended to call at application startup
```

Magnetic field values depend on time and position. Position is specified using WGS-84 geodetic coordinates (latitude, longitude, altitude). Example for calculating magnetic declination:

```java
Field field = IGRFModel.compute(
    40.7128,                    // Latitude (degrees, north positive)
    -74.0060,                   // Longitude (degrees, east positive)
    100.0,                      // Altitude (meters above MSL)
    System.currentTimeMillis()  // Time (epoch milliseconds UTC)
);

double declination = field.declinationDeg;
```

The compute method returns a Field object containing all magnetic field components and derived quantities:

```java
Field field = IGRFModel.compute(latitude, longitude, altitude, epochMillis);

// Access components (all in nanoteslas, except angles in degrees)
field.xNorthNt       // North component
field.yEastNt        // East component  
field.zDownNt        // Down component
field.hHorizontalNt  // Horizontal intensity
field.fTotalNt       // Total field intensity
field.declinationDeg // Declination (magnetic variation)
field.inclinationDeg // Inclination (dip angle)
```

## Coordinates and Output

The library uses WGS-84 geodetic coordinates for input, with latitude from -90° (South) to +90° (North), longitude from -180° to +180° or 0° to 360°, and altitude in meters above the WGS-84 ellipsoid. Output components are provided in the NED (North-East-Down) frame, with X pointing to geographic north, Y pointing east, and Z pointing downward toward Earth's center. All field components are given in nanoteslas (nT), with angles in degrees.

## Time Handling

The library accepts time as epoch milliseconds (UTC) and converts this internally to decimal years for IGRF calculations. You can pass the current time using `System.currentTimeMillis()`, or convert specific dates:

```java
// Specific date
LocalDateTime dt = LocalDateTime.of(2025, 1, 15, 12, 0);
long epochMillis = dt.toInstant(ZoneOffset.UTC).toEpochMilli();
```

For historical dates, the library uses linear interpolation between IGRF epochs. For future dates, secular variation is applied up to 5 years beyond the latest epoch, after which values are clamped with a warning.

## Technical Details

The library implements the International Geomagnetic Reference Field using spherical harmonic expansion to degree and order 13. The magnetic scalar potential $V$ at position $(r, \theta, \phi)$ is:

$$V(r, \theta, \phi) = a \sum_{n=1}^{N} \left(\frac{a}{r}\right)^{n+1} \sum_{m=0}^{n} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] P_n^m(\cos\theta)$$

where $a = 6371.2$ km (reference Earth radius), $g_n^m$ and $h_n^m$ are Schmidt quasi-normalized Gauss coefficients, $P_n^m$ are associated Legendre functions, and $N = 13$ (maximum degree and order). The magnetic field components are derived as $\mathbf{B} = -\nabla V$.

The implementation converts geodetic (WGS-84) to geocentric coordinates, uses linear interpolation between IGRF epochs with secular variation extrapolation up to 5 years beyond the latest epoch, and employs stable recurrence relations for associated Legendre functions. Performance is optimized through thread-local workspace caching and zero-allocation design.

## Logging

The library uses SLF4J for logging. Logging is optional and the library will work silently if no SLF4J binding is present. The library logs:

- **Warnings** when calculations are performed for dates more than 5 years beyond the latest IGRF epoch (coefficients should be updated)
- **Errors** when the coefficient file cannot be found, read, or parsed

To enable logging, add an SLF4J binding to your project like Logback, Log4j2, or slf4j-simple:

```xml
<dependency>
    <groupId>ch.qos.logback</groupId>
    <artifactId>logback-classic</artifactId>
    <version>1.4.11</version>
</dependency>
```

If you copy the source files directly into your project, you'll also need to add the SLF4J API dependency, or remove the logger import and logging statements.

## Coefficients

The library includes embedded IGRF-14 coefficients in the JAR, covering epochs 1900-2030. This is the main benefit of using this library - users don't need to source, manage, or update coefficient files separately. Everything works out of the box. When new IGRF models are released, simply update the library version to get the latest coefficients.

## References

- IGRF Model Documentation: https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
- NOAA Geomagnetic Calculators: https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml
- WGS-84 Specification: https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

This library contains derivative works. See the [NOTICE](NOTICE) file for attribution and original license information for TSAGeoMag.java (public domain, Los Alamos National Laboratory) and GeomagneticField.java (Apache License 2.0, Copyright 2009 The Android Open Source Project).

## Contributing

Contributions are welcome. Please ensure code follows existing conventions, all tests pass, and documentation is updated.

---

**Version 0.0.1-SNAPSHOT** - IGRF-14 model support (epochs 1900-2030)
