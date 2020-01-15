import ggplotnim

proc readSolarModel(fname: string): DataFrame =
  ## first of all parses the given AGSS09 solar model and returns it as a
  ## ggplotnim DataFrame. For efficiency's sake it might be a good idea and
  ## later use a `Table[string, seq[float]]` instead of the dataframe!
  result = toDf(readCsv(fname, sep = ' '))

when isMainModule:
  const solarModel = "./ReadSolarModel/resources/AGSS09_solar_model_stripped.dat"

  let df = readSolarModel(solarModel)
  echo df.pretty(precision = 10)

  # to read a single column, e.g. radius:
  echo df["Radius"]

  # now let's plot radius against temperature colored by density
  ggplot(df, aes("Radius", "Temp", color = "Rho")) +
    geom_line() +
    ggtitle("Radius versus temperature of solar mode, colored by density") +
    ggsave("radius_temp_density.pdf")