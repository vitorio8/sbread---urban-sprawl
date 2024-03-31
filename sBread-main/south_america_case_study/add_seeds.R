ie_historical <-  data.frame(name = "Santa Maria la Antigua del Darien", 
                             year = 1510, lon = -77.021361, lat = 8.21486) %>%
  add_row(name = "Panama City", year = 1519, 
          lon = -79.516667, lat = 8.983333) %>%
  add_row(name = "Cumaná", year = 1521, 
          lon = -64.1675, lat = 10.456389) %>%
  add_row(name = "Santa Marta", year = 1525, 
          lon = -74.205278, lat = 11.241944) %>%
  add_row(name = "Piura", year = 1532, 
          lon = -80.633333, lat = -5.2) %>%
  add_row(name = "Sao Vicente", year = 1532, 
          lon = -46.392222, lat = -23.963333) %>%
  add_row(name = "Cartagena", year = 1533, 
          lon = -75.5, lat = 10.4) %>% 
  add_row(name = "Cuzco", year = 1533, 
          lon = -71.972222, lat = -13.525) %>% 
  add_row(name = "Quito", year = 1534, 
          lon = -78.5125, lat = -0.22) %>% 
  add_row(name = "Trujillo", year = 1534, 
          lon = -79.0288, lat = -8.112) %>%
  add_row(name = "Lima", year = 1535, 
          lon = -77.0375, lat = -12.06) %>%
  add_row(name = "Asunción", year = 1537, 
          lon = -57.633333, lat = -25.3) %>% 
  add_row(name = "Olinda", year = 1537, 
          lon = -34.883333, lat = -8) %>%
  add_row(name = "Bogotá", year = 1538, 
          lon = -74.072222, lat = 4.711111) %>% 
  add_row(name = "Sucre", year = 1538, 
          lon = -65.26, lat = -19.0475) %>%
  add_row(name = "Santiago de Chile", year = 1544, 
          lon = -70.666667, lat = -33.45) %>%
  add_row(name = "Potosí", year = 1545, 
          lon = -65.753333, lat = -19.589167) %>%
  add_row(name = "La Paz", year = 1548, 
          lon = -68.15, lat = -16.5) %>% 
  add_row(name = "Salvador da Bahia", year = 1549, 
          lon = -38.516667, lat = -12.983056) %>%
  add_row(name = "Concepción", year = 1550, 
          lon = -73.051369, lat = -36.828194) %>%
  add_row(name = "Huancavelica", year = 1563, 
          lon = -74.975556, lat = -12.786389) %>% 
  add_row(name = "Rio de Janeiro", year = 1567, 
          lon = -43.205556, lat = -22.911111) %>% 
  add_row(name = "Caracas", year = 1567, 
          lon = -66.903611, lat = 10.480556) %>%
  add_row(name = "Cochabamba", year = 1574, 
          lon =  -66.166667, lat = -17.383333) %>%
  add_row(name = "Buenos Aires", year = 1580, 
          lon = -58.381667, lat = -34.603333) %>% 
  write_csv("data/historical_seeds_ie.csv")