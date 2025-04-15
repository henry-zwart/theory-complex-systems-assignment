#import "@local/uva-game-theory:0.1.0": assignment

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
)

// Display main content (solutions)
#include "solutions.typ"
