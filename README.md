# IGRF Model Library

[![Sonatype Central](https://maven-badges.sml.io/sonatype-central/io.github.tinemuz/igrf-model/badge.png?style=flat)](https://search.maven.org/artifact/io.github.tinemuz/igrf-model)

A lightweight Java library for calculating Earth's magnetic field using the International Geomagnetic Reference Field (IGRF) model for a given location and time.

The library includes embedded IGRF-14 coefficients (epochs 1900-2030). It has no external dependencies (only optional `slf4j-api` for logging) and supports WGS-84 geodetic coordinates.

## Install

### Maven
```xml
<dependency>
    <groupId>io.github.tinemuz</groupId>
    <artifactId>igrf-model</artifactId>
    <version>RELEASE</version>
</dependency>
```

### Gradle
```gradle
implementation 'io.github.tinemuz:igrf-model:+'
```

> **Note:** Check [Maven Central](https://search.maven.org/artifact/io.github.tinemuz/igrf-model) for the latest release version. Snapshot versions are published automatically to Maven Central snapshots repository. For snapshots, add the repository:
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

If you prefer not to add a dependency, you can copy the source files directly into your project. You only need `IGRFModel.java` and the coefficient file `igrfcoeffs.txt` (place it in your resources' folder).

## Example
Simply call `compute()`, coefficients load automatically on first use
```java
Field field = IGRFModel.compute(latitude, longitude, altitude, epochMillis);

// Access components (all in nanoteslas, except angles in degrees)
field.xNorthNt      // North component
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

For historical dates, the library uses linear interpolation between IGRF epochs. For future dates, secular variation is applied up to 5 years beyond the latest epoch, after which values are clamped with a warning to update.

## Technical Details

The library implements the International Geomagnetic Reference Field using spherical harmonic expansion to degree and order 13. The magnetic scalar potential $V$ at position $(r, \theta, \phi)$ is:

$$V(r, \theta, \phi) = a \sum_{n=1}^{N} \left(\frac{a}{r}\right)^{n+1} \sum_{m=0}^{n} \left[ g_n^m \cos(m\phi) + h_n^m \sin(m\phi) \right] P_n^m(\cos\theta)$$

where $a = 6371.2$ km (reference Earth radius), $g_n^m$ and $h_n^m$ are Schmidt quasi-normalized Gauss coefficients, $P_n^m$ are associated Legendre functions, and $N = 13$ (maximum degree and order). The magnetic field components are derived as $\mathbf{B} = -\nabla V$.

## Logging

The library uses SLF4J for logging. Logging is optional and the library will work silently if no SLF4J binding is present.

To enable logging, add an SLF4J binding to your project like Logback, Log4j2, or slf4j-simple:

```xml
<dependency>
    <groupId>ch.qos.logback</groupId>
    <artifactId>logback-classic</artifactId>
    <version>LATEST</version>
</dependency>
```

## References

- IGRF Model Documentation: https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
- NOAA Geomagnetic Calculators: https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml
- WGS-84 Specification: https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

This library contains derivative works. See the [NOTICE](NOTICE) file for attribution and original license information for TSAGeoMag.java (public domain, Los Alamos National Laboratory) and GeomagneticField.java (Apache License 2.0, Copyright 2009 The Android Open Source Project).

## Contributing

Contributions are welcome. Please ensure tests pass, and documentation is updated.

---

**IGRF-14 model support** (epochs 1900-2030)
