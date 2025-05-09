#import "../uva-game-theory/src/lib.typ": assignment 

#set text(lang: "EN", region: "GB")

#let authors = (
  (name: "Henry Zwart", id: "15393879"),
)

#let course = (
  lecturer: "Dr. Clélia de Mulatier",
  name: "Theory of Complex Systems",
  code: "5284THCS6Y",
)

#show: assignment.with(
  title: "Assignment",
  course: course,
  authors: authors,
  font: "New Computer Modern",
  fontsize: 9pt,
  margin_size: 2cm,
)

// Display main content (solutions)
#include "solutions.typ"
