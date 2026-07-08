// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// https://astro.build/config
export default defineConfig({
	site: 'https://dbretina.github.io',
	base: '/DBRetina/',
	devToolbar: { enabled: false },
	integrations: [
		starlight({
			title: 'DBRetina',
			description:
				'Pairwise similarity between supergroups (diseases, drugs, pathways) from their shared features.',
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/DBRetina/DBRetina' },
			],
			lastUpdated: true,
			sidebar: [
				{
					label: 'Getting started',
					items: [
						{ label: 'Installation', slug: 'start/installation' },
						{ label: 'Quickstart', slug: 'start/quickstart' },
						{ label: 'Core concepts', slug: 'start/concepts' },
					],
				},
				{
					label: 'Commands',
					items: [
						{
							label: 'Build & index',
							items: [
								{ label: 'index', slug: 'commands' },
								{ label: 'merge', slug: 'commands/merge' },
							],
						},
						{
							label: 'Core analysis',
							items: [
								{ label: 'pairwise', slug: 'commands/pairwise' },
								{ label: 'query', slug: 'commands/query' },
								{ label: 'cluster', slug: 'commands/cluster' },
								{ label: 'export', slug: 'commands/export' },
							],
						},
						{
							label: 'Gene-set analysis',
							items: [
								{ label: 'dedup', slug: 'commands/dedup' },
								{ label: 'modularity', slug: 'commands/modularity' },
								{ label: 'setcov', slug: 'commands/setcov' },
							],
						},
						{
							label: 'Graph & network',
							items: [
								{ label: 'bipartite', slug: 'commands/bipartite' },
								{ label: 'graph', slug: 'commands/graph' },
							],
						},
						{
							label: 'Explore an index',
							items: [
								{ label: 'geneinfo', slug: 'commands/geneinfo' },
								{ label: 'Quick lookups', slug: 'commands/lookups' },
							],
						},
						{
							label: 'Network analysis',
							items: [
								{ label: 'connect', slug: 'commands/connect' },
								{ label: 'module', slug: 'commands/module' },
								{ label: 'enrich', slug: 'commands/enrich' },
							],
						},
					],
				},
				{
					label: 'Examples',
					items: [
						{ label: 'Overview', slug: 'examples' },
						{ label: 'Shared autoimmune core', slug: 'examples/autoimmune-diseases' },
						{ label: 'Comorbidity, not cause', slug: 'examples/autism-comorbidity' },
						{ label: 'One name, two mechanisms', slug: 'examples/diabetes-types' },
						{ label: "Alzheimer's and set size", slug: 'examples/alzheimer-metric' },
						{ label: 'Pathway database concordance', slug: 'examples/database-concordance' },
						{ label: 'Hallmark independence', slug: 'examples/hallmarks' },
						{ label: 'Drug mechanism from targets', slug: 'examples/drug-repurposing' },
						{ label: 'Disease to drug, five databases', slug: 'examples/multidb-rheumatoid' },
						{ label: "Recovering cancer's machinery", slug: 'examples/enrich-cancer-hallmarks' },
						{ label: 'Cancer across five databases', slug: 'examples/enrich-crossdb-cancer' },
					],
				},
				{
					label: 'Reference',
					items: [
						{ label: 'File formats', slug: 'reference/formats' },
						{ label: 'Metrics', slug: 'reference/metrics' },
						{ label: 'Command chaining', slug: 'reference/chaining' },
					],
				},
			],
		}),
	],
});
