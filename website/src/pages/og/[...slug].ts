import { getCollection } from 'astro:content';
import { OGImageRoute } from 'astro-og-canvas';

// One OpenGraph card per docs page, rendered at build time on the brand gradient.
const entries = await getCollection('docs');
const pages = Object.fromEntries(entries.map((entry) => [entry.id || 'index', entry]));

export const { getStaticPaths, GET } = await OGImageRoute({
	param: 'slug',
	pages,
	getImageOptions: (_path, entry) => ({
		title: entry.data.title,
		description:
			entry.data.description ??
			'Similarity networks for diseases, drugs, and pathways, from the features they share.',
		bgGradient: [
			[20, 16, 42],
			[16, 40, 44],
		],
		border: { color: [124, 92, 255], width: 16, side: 'inline-start' },
		padding: 70,
		font: {
			title: { color: [255, 255, 255], size: 74, weight: 'Bold', lineHeight: 1.1 },
			description: { color: [190, 198, 214], size: 34, lineHeight: 1.4 },
		},
	}),
});
