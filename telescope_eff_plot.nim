import ggplotnim

let fname2 = "extracted_from_aznar2015_llnl_telescope_eff_plot.csv"
let dfEnergyEff = toDf(readCsv(fname2, sep = ','))
  .mutate(fn {"Energies [keV]" ~ `xVals` * 8.0 + 1.0}, 
          fn {"Effective Area [cm^2]" ~ `yVals` * 8.5})

ggplot(dfEnergyEff, aes("Energies [keV]", "Effective Area [cm^2]")) +
    geom_line() +
    ggtitle("The telescope energy efficiency") +
    ggsave("EnergyEff.pdf")
echo dfEnergyEff