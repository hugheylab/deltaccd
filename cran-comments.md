## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

  > checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Jake Hughey <jakejhughey@gmail.com>'
    
    New submission
    
    Possibly misspelled words in DESCRIPTION:
      Shilts (10:32)
      al (10:42)
      et (10:39)
      transcriptome (8:56)
  
  > checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
      'lastMiKTeXException'
  
  0 errors ✓ | 0 warnings ✓ | 2 notes x

Notes:
  - This is a resubmission of the first time this package was sent to CRAN, and the words are names and correctly spelled terms, so ignore this.
  - The second note only occurs on the Windows Server rhub environment, and from what I have seen about these types of notes they do not occur when building and checking on CRAN.

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/deltaccd/actions).

## Resubmission notes
Email response to initial submission:
"Thanks,

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.

If you use a package which is only needed in examples, please list it in
'Suggests' and wrap these examples in if(requireNamespace("pkgname")){}
instead.


Please fix and resubmit.

Best,
Gregor Seyer"

  - We updated the examples and documentation based on this feedback.

## Downstream dependencies
There are no downstream dependencies for deltaccd.

## Tests
When checking locally, tests take about 2.1s total. If this is too long for tests let us know and we can skip some tests on CRAN.