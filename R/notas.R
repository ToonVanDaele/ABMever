
# Nog te doen:

# - 1 januari is ok als startdatum?  =>  1 april start  (generiek implemeteren)
# -> bekijken of de functie voor de initiële populatie aangepast moet worden!!

# set_init_pop   functie aan te passen aan Sm in de vorm van een matrix (nu enkel
# als vector length 3)

# maand (kalendermaand) mee bewaren bij afschot en aantallen

# abm simulaties -> overwegen om leeftijdsklassen volledig te laten vallen
# en enkel te werken met de effectieve leeftijd (met een absoluut maximum vb. 15 jaar)
# Dit kan het selectieproces en de indexen sterk vereenvoudigen
# Voorbereidende verwerking moet dan voor elke parameter een vector maken
# met voor elke maand een waarde voor overleving, fertiliteit, ...

# afschot tweestaps relatief afschot. bv. 0.1 van totale populatie. Telkens
# per klasse volgens vaste verhouding.

# is numboar na of voor jacht in het opgegeven tijdstap?
# volgorde: reproduce -> mortality -> hunting -> (aging)
# bekijken of er effect is  -> waarom vershil bij scenario example.

# df_harvest -> volledige populatie informatie (boar - who dies) bijhouden.
# + sim, scenario en tijdstap (tijd + kalendermaand) + leeftijd + geslacht
# zowel populatie als harvest volledige populatie.



# Overige nota's

# - voorlopig 3 leeftijdsklassen voor populatie behouden (niet naar 4 categoriën)
# -> ? generiek implementeren

# -> figuur maken van maximale leeftijd met en zonder jacht -> maximum ongeveer 15 jaar?


# Een afschot van een bepaald jaar -> hoe hetzelfde afschot verdelen in het jaar
# ttz een jachtseizoen introduceren.
#
# Heeft het zin om alle afschot voor de geboortepiek te concentreren?
#
# zomerafschot is om schade te vermijden, minder om de populatie te reguleren
# ? wat is de impact van het zomerafschot op

# Waar je controle over hebt (qua jacht) zijn de absolute aantallen,
# de relatieve verdeling over het jaar en tussen de leeftijdsklassen.
# Er is geen controle over de ratios afschot t.a.v. de werkelijke populatie.

# Beginnen we met de stable stage zonder jacht of stable stage met jacht?

# Er is nu geen maximale leeftijd voor overleving. ? maximum instellen?
# Of overleving leeftijdsafhankelijk? Briederman?? -> nog te bekijken (Jim)
# repercussie ->

