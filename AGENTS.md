# Agent workspace rules

- Keep all agent-created files inside this repository. Do not write code, data, or other task artifacts to `/tmp` or another external directory.
- Put exploratory code, one-off test programs, generated data, reports, logs, and temporary build output in `.agentic_artifacts/`. Create it if it is absent. Its contents are ignored by Git.
- Put requested, lasting project changes in the normal source and documentation directories so they can be reviewed and committed.
- Do not move, delete, or overwrite existing project files as part of workspace cleanup unless the user asks.
