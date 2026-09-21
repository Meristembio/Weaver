# Protocols

The `protocols` app manages reusable laboratory preparation data. It is separate from the DNA inventory but can be reached from the main navigation through `Recipes` and `Components`.

## Objects and Ownership

- A reactive is a solid or liquid reagent owned by a user. Record its name, state, molar mass for solids when needed, concentration and unit for liquids, and whether it is autoclavable.
- A component connects one owned reactive to a target concentration and unit.
- A recipe has a name, four-character code, optional pH, description, categories, required components, optional components, and projects to which it is shared.
- A variant belongs to a recipe and has its own name, four-character code, pH, description, and optional components chosen from that recipe.

Components can use M, mM, mg/ml, volume/volume percent, or weight/volume percent units. A recipe code and variant code must be four characters when supplied, and recipe codes are unique.

## Create the Recipe Data

1. Open `Components` and create the reactive records you own. For a solid, provide molar mass if later calculations will use molar units. For a liquid, record the stock concentration and unit. Mark non-autoclavable reagents so the preparation result can warn about them.
2. Create a component and select one of your reactives, then enter its target concentration and unit. A component is not available to another user unless that user owns it or it is otherwise exposed by the current deployment.
3. Open `Recipes > Create`, enter the recipe name and optional code, pH, description, categories, required components, optional components, and projects in `Shared to project`. Only projects where you have write or admin access are available for sharing.
4. Save the recipe and open its detail page. Create a variant from the recipe when a preparation needs a defined alternative; its optional-component list is limited to components already marked optional on the parent recipe.

## Calculate a Preparation

Open a recipe or variant and enter a quantity, unit (`ml` or `lt`), and optional concentration multiplier in the preparation form. Weaver calculates the required mass or volume for each component, converts small quantities to mg, ml, or ul where appropriate, and lists a warning when it cannot calculate a value. Typical errors include missing molar mass, incompatible mass/volume conversions, a stock concentration that is equal to or lower than the target, or an attempt to convert a solid into volume/volume units. The result also indicates when one or more reactives are not autoclavable.

For a recipe, the calculation uses required components. For a variant, it uses the parent recipe's required components plus the variant's selected optional components. Check every calculated amount and every warning before preparing the solution.

## Permissions and Labels

Recipe editing is owner-only. A recipe can be viewed when it is owned by the user, shared to one of the user's projects, or visible through the current project's access rules. Components and reactives are owned records, so create or edit them under the account that should maintain the data. Recipe and variant label pages are read views derived from the saved records; they do not replace the underlying recipe or variant.
