# swiftc fails on the default SDK, so the Vision OCR witness silently produced empty files

- Date: 2026-10-11
- Agent: corpus-henon-opus
- Seen before: no (the 1974 Hénon transcription built the same vision_ocr.swift earlier, before the SDK 27 install)

What happened: `swiftc vision_ocr.swift` (the Vision witness script filed beside the Hénon 1974 PDF) failed with "this SDK is not supported by the compiler" (default SDK MacOSX27.0, compiler Swift 6.3.3). My loop ran `vision_ocr` without checking that it had been built, so all 23 witness files were empty and the first count showed Vision = 0 for every number.
Workaround: `swiftc -O -sdk /Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk vision_ocr.swift -o vision_ocr` builds and runs.
Suggested fix: pass the -sdk flag in any witness recipe, and make the witness loop fail on a missing binary or an empty output file (an all-zero witness is an absent instrument, not a negative).

Disposition (2026-10-11 19:58, lead): Promoted (coordination.md numerics rules: a witness loop fails on a missing binary or an empty output file; an all-zero witness is an absent instrument) and Known environment (default SDK 27 vs swiftc; build with -sdk MacOSX26.5.sdk). The machine note line is proposed to the owner (user-level instructions).
