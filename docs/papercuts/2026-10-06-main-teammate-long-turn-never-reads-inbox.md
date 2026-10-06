# A teammate in one long turn never reads the lead's replies, even to its own questions

- Date: 2026-10-06
- Agent: main
- Seen before: yes (processed/2026-10-05-main-team-lead-address-unread.md, the first Fable reviewer's 5 unread messages)

What happened: twobody-gen-opus worked in one continuous turn for hours, sending "ruling needed" questions to main. All 10 lead replies since about 13:14 AEDT sat unread in its inbox file (read=false), so it kept re-asking and proceeding on its own recommendations.
Workaround: none needed for safety (it proceeded only on safe steps); lead noticed via the inbox file.
Suggested fix: briefs say "after any message that asks the lead for a ruling, END YOUR TURN so the reply can arrive; resume when it does". Lead checks the teammate's inbox file for read=false when a teammate repeats a question.
