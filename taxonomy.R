# Taxonomy - create lookup tables

# Load taxonomy
taxonomy <- load_taxa(server=getOption("fishbase"))

# Create species lookup table
species_lookup <- taxonomy %>%
 select(Species, SpecCode, Order, Class, SuperClass) %>%
 rename(accepted_name = Species) %>%
 distinct()

# Create genus lookup table
genus_lookup <- taxonomy %>%
  select(Genus, SpecCode, Order, Class, SuperClass) %>%
  rename(accepted_name = Genus) %>%
  distinct()

# Create family lookup table
# family_lookup <- taxonomy %>%
# select(Family, SpecCode, Order, Class, SuperClass) %>%
# rename(accepted_name = Family) %>%
# distinct()

# Combine
taxonomy_lookup <- bind_rows(species_lookup, genus_lookup) %>%
 distinct(accepted_name, .keep_all = TRUE) %>%
 rename(fishbase_code = SpecCode, order = Order, class = Class, superclass = SuperClass)

# Check each name only occurs once
taxonomy_lookup %>%
 count(accepted_name) %>%
 filter(n > 1)

# Save as CSV
write_csv(taxonomy_lookup, "taxonomy_lookup.csv")


