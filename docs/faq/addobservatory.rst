.. doctest-skip-all

.. _addobservatory:

********************************************************
How can I add observatories to the observatory database?
********************************************************

`astroplan` has a convenience methods for creating an
`~astropy.coordinates.EarthLocation` object at the location of an observatory
– see `~astroplan.get_site`, `~astroplan.get_site_names` and
`~astroplan.add_site` for more on what you can do with the observatory
database).

The observatory database is stored as a JSON file in the `astroplan/data`
directory. If you would like to add an observatory to the observatory database
for distribution in the `astroplan` source code, use the
`~astroplan.new_site_info_to_json` function to output the observatory's
information into a string with the proper JSON format. Then you can copy and
paste that JSON string with the new observatory's information into the
astroplan/data/observatories.json` file and submit a
`pull request on GitHub <https://github.com/astropy/astroplan/pulls>`_ to add
that observatory to the source code.
