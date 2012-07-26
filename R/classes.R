setClass("experiment",
         representation=representation(
           result.file="character",
           proteins="list",
           parameters="rTParam",
           unused.parameters="vector",
           xtandem.version="character",
           start.time="character",
           source.paths="vector",
           source.paths.descriptions="vector",
           estimated.false.positive="integer",
           quality.values="vector",
           total.peptides.used="integer",
           total.proteins.used="integer",
           total.spectra.assigned="integer",
           input.spectra="integer",
           input.models="integer",

           partial.cleavages="integer",
           point.mutation="integer",
           unanticipated.cleavages="integer",
           potential.C.terminii="integer",
           potential.N.terminii="integer",

           initial.modelling.time="numeric",
           initial.spectrum.modelling.time="numeric",
           load.sequence.models.time="numeric",
           refinement.per.spectrum.time="numeric"
           )
         )

setClass("protein",
         representation=representation(
           expect.value="numeric",
           file="character",
           label="character",
           sequence="character",
           uid="integer",
           peptides="list"
           )
         )

setClass("peptides",
         representation=representation(
           expect.value="numeric",
           tandem.score="numeric",
           mh="numeric",
           delta="numeric",
           peak.count="integer",
           missed.cleavage="integer",
           start.position="integer",
           end-position="integer",
           sequence="character",
           spectrum="spectrum",
           ptm="list"
           )
         )

setClass("ptm",
         representation=representation(
           type="character",
           position="integer",
           mass.change="numeric"
           )
         )

setClass("spectrum",
         representation=representation(
           id="integer",
           mh="numeric",
           sumI="numeric",
           maxI="numeric",
           fI="numeric"
           )
         )
          
          
           
           
          
